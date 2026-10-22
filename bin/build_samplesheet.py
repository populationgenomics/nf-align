#!/usr/bin/env python3

"""
Build the alignment sample sheet for one Metamist project and sequencing type from a list of cohorts.

Run by the BUILD_SAMPLESHEET process at the start of main.nf (see DESIGN.md). The input rules mirror
cpg-flow (cpg_flow.inputs, cpg_flow.metamist), restricted to a single project:
- cohorts are looked up within the project; missing, inactive and empty cohorts fail the run
- sequencing groups are de-duplicated across cohorts, and must belong to the project
- SGs that aren't short-read, or aren't of the requested sequencing type, are skipped with a warning
- reads are combined across an SG's assays, each assay contributing one FASTQ pair (or the CRAM/BAM)
- an SG with a completed, active 'cram' analysis is already aligned, and gets one 'aligned' row
- read locations are used exactly as recorded in Metamist

    build_samplesheet.py --project <dataset>[-test] --sequencing-type genome \\
        --cohorts COH0001 COH0002 --output samplesheet.tsv
"""

import csv
import logging
from argparse import ArgumentParser
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'kind', 'n_units', 'unit', 'r1', 'r2', 'reads', 'reads_index', 'reads_ref']
SEQUENCING_TYPES = {'genome', 'exome'}

# one round trip: cohort status, SGs, assay metadata, and any existing CRAM analyses
COHORTS_QUERY = """
query BuildSamplesheet($project: String!, $cohorts: [String!]!) {
    project(name: $project) {
        cohorts(id: {in_: $cohorts}) {
            id
            status
            sequencingGroups {
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
        }
    }
}
"""


class SamplesheetError(Exception):
    """Invalid cohort or sequencing group metadata."""


def fetch_cohorts(project: str, cohort_ids: list[str]) -> list[dict]:
    """Query the project's cohorts, failing on any that are missing (from this project), inactive, or empty."""
    cohorts = query(gql(COHORTS_QUERY), variables={'project': project, 'cohorts': cohort_ids})['project']['cohorts']

    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts in {project} - missing: {missing}, inactive: {inactive}, empty: {empty}',
        )
    return cohorts


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
    ref = (meta.get('reference_assembly') or {}).get('location', '') if reads_type == 'cram' else ''
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


def build_rows(cohorts: list[dict], project: str, sequencing_type: str) -> list[dict]:
    """Rows for every in-scope SG. Reports every bad SG at once, not one per run."""
    rows: list[dict] = []
    seen: set[str] = set()
    errors: list[str] = []
    for cohort in cohorts:
        for sg in cohort['sequencingGroups']:
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
        raise SamplesheetError(f'No short-read {sequencing_type} sequencing groups in these cohorts')
    return rows


def write_sheet(rows: list[dict], path: str) -> None:
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', restval='')
        writer.writeheader()
        writer.writerows(rows)


def main(project: str, sequencing_type: str, cohort_ids: list[str], output: str) -> None:
    rows = build_rows(fetch_cohorts(project, sorted(set(cohort_ids))), project, sequencing_type)
    write_sheet(rows, output)
    n_sgs = len({row['sg_id'] for row in rows})
    n_aligned = sum(row['kind'] == 'aligned' for row in rows)
    logger.info(f'{project} {sequencing_type}: {n_sgs} SGs ({n_aligned} already aligned), {len(rows)} rows -> {output}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, help='Metamist project: the dataset, or <dataset>-test')
    parser.add_argument('--sequencing-type', required=True, choices=sorted(SEQUENCING_TYPES))
    parser.add_argument('--cohorts', nargs='+', required=True, help='Metamist cohort IDs')
    parser.add_argument('--output', required=True, help='Sample sheet TSV to write')
    args = parser.parse_args()
    main(args.project, args.sequencing_type, args.cohorts, args.output)
