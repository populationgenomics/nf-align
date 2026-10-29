#!/usr/bin/env python3

"""
Build the alignment sample sheet for one Metamist project and sequencing type, from a list of cohorts or a
list of sequencing group IDs.

Run by the BUILD_SAMPLESHEET process at the start of main.nf (see DESIGN.md). The input rules mirror
cpg-flow (cpg_flow.inputs, cpg_flow.metamist), restricted to a single project:
- cohorts are looked up within the project; missing, inactive and empty cohorts fail the run
- sequencing groups are de-duplicated across cohorts, and must belong to the project
- SGs from cohorts that aren't short-read, or aren't of the requested sequencing type, are skipped with a warning
- SG IDs are looked up within the project; missing SGs, and SGs that aren't short-read of the requested
  sequencing type, fail the run, as they were asked for by name
- reads are combined across an SG's assays, each assay contributing one FASTQ pair (or the CRAM/BAM)
- an SG with a completed, active 'cram' analysis is already aligned, and gets one 'aligned' row
- read locations are used exactly as recorded in Metamist

    build_samplesheet.py --project <dataset>[-test] --sequencing-type genome \\
        --cohorts COH0001 COH0002 --output samplesheet.tsv
    build_samplesheet.py --project <dataset>[-test] --sequencing-type genome \\
        --sgids CPG0001 CPG0002 --output samplesheet.tsv
"""

import csv
import logging
from argparse import ArgumentParser
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'kind', 'n_units', 'unit', 'r1', 'r2', 'reads', 'reads_index', 'reads_ref']
SEQUENCING_TYPES = {'genome', 'exome'}

# the SG fields both queries fetch: assay metadata, and any existing CRAM analyses
SG_FIELDS = """
fragment SgFields on GraphQLSequencingGroup {
    id
    type
    technology
    sample {
        project {
            name
        }
    }
    assays {
        meta
    }
    analyses(type: {eq: "cram"}, status: {eq: COMPLETED}, active: {eq: true}) {
        id
        output
    }
}
"""

# one round trip: cohort status, and each cohort's SGs
COHORTS_QUERY = (
    """
query BuildSamplesheetFromCohorts($project: String!, $cohorts: [String!]!) {
    project(name: $project) {
        cohorts(id: {in_: $cohorts}) {
            id
            status
            sequencingGroups {
                ...SgFields
            }
        }
    }
}
"""
    + SG_FIELDS
)

# one round trip: the requested SGs, queried directly
SGS_QUERY = (
    """
query BuildSamplesheetFromSgs($project: String!, $sgids: [String!]!) {
    project(name: $project) {
        sequencingGroups(id: {in_: $sgids}) {
            ...SgFields
        }
    }
}
"""
    + SG_FIELDS
)


class SamplesheetError(Exception):
    """Invalid cohorts, sequencing group IDs, or sequencing group metadata."""


def fetch_cohort_sgs(project: str, cohort_ids: list[str]) -> list[dict]:
    """
    The SGs of the project's cohorts, in cohort order and possibly repeated. Fails on any cohort that is
    missing (from this project), inactive, or empty.
    """
    cohorts = query(gql(COHORTS_QUERY), variables={'project': project, 'cohorts': cohort_ids})['project']['cohorts']

    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts in {project} - missing: {missing}, inactive: {inactive}, empty: {empty}',
        )
    return [sg for cohort in cohorts for sg in cohort['sequencingGroups']]


def fetch_sgs(project: str, sg_ids: list[str], sequencing_type: str) -> list[dict]:
    """
    The project's SGs with these IDs. Each was asked for by name, so one that is missing (from this project),
    or isn't short-read of this sequencing type, fails the run rather than being skipped.
    """
    sgs = query(gql(SGS_QUERY), variables={'project': project, 'sgids': sg_ids})['project']['sequencingGroups']

    missing = sorted(set(sg_ids) - {sg['id'] for sg in sgs})
    wrong_type = sorted(
        f'{sg["id"]} ({sg["technology"]} {sg["type"]})'
        for sg in sgs
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type
    )
    if missing or wrong_type:
        raise SamplesheetError(
            f'Unusable SGs in {project} - missing: {missing}, not short-read {sequencing_type}: {wrong_type}',
        )
    return sgs


def combine_assay_meta(assays: list[dict]) -> dict:
    """Each assay's 'reads' becomes one entry in a combined list; other keys are last-wins."""
    combined: dict = {}
    for assay in assays:
        for key, value in assay['meta'].items():
            if key == 'reads':
                combined.setdefault('reads', []).append(value)
            else:
                combined[key] = value
    return combined


def alignment_file_rows(sg_id: str, meta: dict) -> list[dict]:
    """One row for a single CRAM/BAM input, with its index (if any) and the CRAM's reference."""
    reads_type, reads = meta['reads_type'], meta['reads']
    if len(reads) != 1:
        raise SamplesheetError(f'{sg_id}: expected one {reads_type} input, got {len(reads)}')
    location = reads[0]['location']
    if not location.endswith(f'.{reads_type}'):
        raise SamplesheetError(f'{sg_id}: reads_type is {reads_type}, but the file is {location}')
    index = ''
    if secondary := reads[0].get('secondaryFiles'):
        index = secondary[0]['location']
        if not index.endswith('.crai' if reads_type == 'cram' else '.bai'):
            raise SamplesheetError(f'{sg_id}: unexpected index {index} for {location}')
    ref = ''
    if reads_type == 'cram':
        # a CRAM decodes only against the reference it was written with, so there's no falling back to ref_fasta
        ref = (meta.get('reference_assembly') or {}).get('location') or ''
        if not ref:
            raise SamplesheetError(f'{sg_id}: "reference_assembly.location" is required for cram assays')
    return [
        {
            'sg_id': sg_id,
            'kind': reads_type,
            'n_units': 1,
            'unit': 1,
            'reads': location,
            'reads_index': index,
            'reads_ref': ref,
        },
    ]


def fastq_rows(sg_id: str, meta: dict) -> list[dict]:
    """One row per FASTQ / FASTQ.ora pair; ORA rows carry the decompression reference."""
    reads_type, reads = meta['reads_type'], meta['reads']
    ora_ref = ''
    if reads_type == 'fastq_ora':
        ora_ref = (meta.get('ora_reference') or {}).get('location') or ''
        if not ora_ref:
            raise SamplesheetError(f'{sg_id}: "ora_reference.location" is required for fastq_ora assays')
    rows = []
    for unit, pair in enumerate(reads, start=1):
        if len(pair) != 2:  # noqa: PLR2004
            raise SamplesheetError(f'{sg_id}: expected an R1/R2 pair per assay, got {len(pair)} files')
        rows.append(
            {
                'sg_id': sg_id,
                'kind': reads_type,
                'n_units': len(reads),
                'unit': unit,
                'r1': pair[0]['location'],
                'r2': pair[1]['location'],
                'reads_ref': ora_ref,
            },
        )
    return rows


def rows_for(sg: dict) -> list[dict]:
    """Sample sheet rows for one SG: one 'aligned' row if it has a CRAM, else one per FASTQ lane, else one."""
    sg_id = sg['id']
    # the latest CRAM analysis wins
    if crams := sorted((a for a in sg['analyses'] if a.get('output')), key=lambda a: a['id']):
        cram = crams[-1]['output']
        return [
            {'sg_id': sg_id, 'kind': 'aligned', 'n_units': 1, 'unit': 1, 'reads': cram, 'reads_index': f'{cram}.crai'},
        ]

    meta = combine_assay_meta(sg['assays'])
    if not meta.get('reads'):
        raise SamplesheetError(f'{sg_id}: no CRAM analysis, and no "reads" in assay meta')
    if meta.get('reads_type') in {'cram', 'bam'}:
        return alignment_file_rows(sg_id, meta)
    if meta.get('reads_type') in {'fastq', 'fastq_ora'}:
        return fastq_rows(sg_id, meta)
    raise SamplesheetError(f'{sg_id}: unsupported reads_type {meta.get("reads_type")!r}')


def build_rows(sgs: list[dict], project: str, sequencing_type: str) -> list[dict]:
    """Rows for every in-scope SG, each SG once. Reports every bad SG at once, not one per run."""
    rows: list[dict] = []
    seen: set[str] = set()
    errors: list[str] = []
    for sg in sgs:
        if sg['id'] in seen:
            continue
        seen.add(sg['id'])
        if (sg_project := sg['sample']['project']['name']) != project:
            errors.append(f'{sg["id"]}: belongs to project {sg_project}, not {project}')
            continue
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type:
            logger.warning(
                f'{sg["id"]}: skipping {sg["technology"]} {sg["type"]}, this run is short-read {sequencing_type}'
            )
            continue
        try:
            rows.extend(rows_for(sg))
        except SamplesheetError as error:
            errors.append(str(error))
    if errors:
        raise SamplesheetError('\n'.join(errors))
    if not rows:
        raise SamplesheetError(f'No short-read {sequencing_type} sequencing groups to align')
    return rows


def write_sheet(rows: list[dict], path: str) -> None:
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', restval='')
        writer.writeheader()
        writer.writerows(rows)


def main(
    project: str,
    sequencing_type: str,
    cohort_ids: list[str] | None,
    sg_ids: list[str] | None,
    output: str,
) -> None:
    if cohort_ids:
        sgs = fetch_cohort_sgs(project, sorted(set(cohort_ids)))
    else:
        sgs = fetch_sgs(project, sorted(set(sg_ids or [])), sequencing_type)
    rows = build_rows(sgs, project, sequencing_type)
    write_sheet(rows, output)
    n_sgs = len({row['sg_id'] for row in rows})
    n_aligned = sum(row['kind'] == 'aligned' for row in rows)
    logger.info(f'{project} {sequencing_type}: {n_sgs} SGs ({n_aligned} already aligned), {len(rows)} rows -> {output}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, help='Metamist project: the dataset, or <dataset>-test')
    parser.add_argument('--sequencing-type', required=True, choices=sorted(SEQUENCING_TYPES))
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--cohorts', nargs='+', help='Metamist cohort IDs')
    inputs.add_argument('--sgids', nargs='+', help='Metamist sequencing group IDs, queried directly')
    parser.add_argument('--output', required=True, help='Sample sheet TSV to write')
    args = parser.parse_args()
    main(args.project, args.sequencing_type, args.cohorts, args.sgids, args.output)
