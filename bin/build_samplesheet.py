#!/usr/bin/env python3

"""
Build the alignment sample sheet for one sequencing type, from a list of cohorts or a list of sequencing
group IDs. The SGs can come from any number of Metamist projects, and each row records its SG's dataset.

Run by the BUILD_SAMPLESHEET process at the start of main.nf (see DESIGN.md). The input rules mirror
cpg-flow (cpg_flow.inputs, cpg_flow.metamist):
- cohorts are looked up across projects; missing (or not visible to this run), inactive and empty cohorts
  fail the run
- sequencing groups are de-duplicated across cohorts
- an SG's dataset is its sample's project, without any -test suffix. Its project must match the access
  level: a -test project at access level test, and any other project otherwise. Mismatches fail the run
- SGs from cohorts that aren't short-read, or aren't of the requested sequencing type, are skipped with a warning
- SG IDs are looked up across projects; missing SGs (or ones not visible to this run), and SGs that aren't
  short-read of the requested sequencing type, fail the run, as they were asked for by name
- reads are combined across an SG's assays, each assay contributing one FASTQ pair (or the CRAM/BAM)
- an SG with a completed, active 'cram' analysis is already aligned, and gets one 'aligned' row. That row
  carries the CRAM's reference (if its analysis records one), and its 'somalier' fingerprint (if it has a
  completed, active one). main.nf fingerprints an aligned SG that has none
- read locations are used exactly as recorded in Metamist

    build_samplesheet.py --access-level standard --sequencing-type genome \\
        --cohorts COH0001 COH0002 --output samplesheet.tsv
    build_samplesheet.py --access-level test --sequencing-type genome \\
        --sgids CPG0001 CPG0002 --output samplesheet.tsv
"""

import csv
import logging
from argparse import ArgumentParser
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'dataset', 'kind', 'n_units', 'unit', 'r1', 'r2', 'reads', 'reads_index', 'reads_ref', 'somalier']
SEQUENCING_TYPES = {'genome', 'exome'}
ACCESS_LEVELS = {'full', 'standard', 'test'}

# the SG fields both queries fetch: assay metadata, and any existing CRAM and somalier analyses
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
    analyses(type: {in_: ["cram", "somalier"]}, status: {eq: COMPLETED}, active: {eq: true}) {
        id
        type
        output
        meta
    }
}
"""

# one round trip: cohort status, and each cohort's SGs. Top-level, so the cohorts can be in any project
COHORTS_QUERY = (
    """
query BuildSamplesheetFromCohorts($cohorts: [String!]!) {
    cohorts(id: {in_: $cohorts}) {
        id
        status
        sequencingGroups {
            ...SgFields
        }
    }
}
"""
    + SG_FIELDS
)

# one round trip: the requested SGs, queried directly, in any project
SGS_QUERY = (
    """
query BuildSamplesheetFromSgs($sgids: [String!]!) {
    sequencingGroups(id: {in_: $sgids}) {
        ...SgFields
    }
}
"""
    + SG_FIELDS
)


class SamplesheetError(Exception):
    """Invalid cohorts, sequencing group IDs, or sequencing group metadata."""


def fetch_cohort_sgs(cohort_ids: list[str]) -> list[dict]:
    """
    The SGs of these cohorts, in cohort order and possibly repeated. Fails on any cohort that is missing
    (or not visible to this run's credentials), inactive, or empty.
    """
    cohorts = query(gql(COHORTS_QUERY), variables={'cohorts': cohort_ids})['cohorts']

    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts - missing or not visible to this run: {missing}, inactive: {inactive}, empty: {empty}',
        )
    return [sg for cohort in cohorts for sg in cohort['sequencingGroups']]


def fetch_sgs(sg_ids: list[str], sequencing_type: str) -> list[dict]:
    """
    The SGs with these IDs. Each was asked for by name, so one that is missing (or not visible to this run's
    credentials), or isn't short-read of this sequencing type, fails the run rather than being skipped.
    """
    sgs = query(gql(SGS_QUERY), variables={'sgids': sg_ids})['sequencingGroups']

    missing = sorted(set(sg_ids) - {sg['id'] for sg in sgs})
    wrong_type = sorted(
        f'{sg["id"]} ({sg["technology"]} {sg["type"]})'
        for sg in sgs
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type
    )
    if missing or wrong_type:
        raise SamplesheetError(
            f'Unusable SGs - missing or not visible to this run: {missing}, '
            f'not short-read {sequencing_type}: {wrong_type}',
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


def latest_analysis(analyses: list[dict], analysis_type: str) -> dict | None:
    """The latest analysis of this type that has an output, if any."""
    found = sorted((a for a in analyses if a['type'] == analysis_type and a.get('output')), key=lambda a: a['id'])
    return found[-1] if found else None


def rows_for(sg: dict) -> list[dict]:
    """Sample sheet rows for one SG: one 'aligned' row if it has a CRAM, else one per FASTQ lane, else one."""
    sg_id = sg['id']
    if cram := latest_analysis(sg['analyses'], 'cram'):
        somalier = latest_analysis(sg['analyses'], 'somalier')
        return [
            {
                'sg_id': sg_id,
                'kind': 'aligned',
                'n_units': 1,
                'unit': 1,
                'reads': cram['output'],
                'reads_index': f'{cram["output"]}.crai',
                # recorded by nf-align's registration; older CRAMs leave main.nf to assume --aligned_ref_fasta
                'reads_ref': (cram.get('meta') or {}).get('ref_fasta', ''),
                # empty: main.nf fingerprints the CRAM
                'somalier': somalier['output'] if somalier else '',
            },
        ]

    meta = combine_assay_meta(sg['assays'])
    if not meta.get('reads'):
        raise SamplesheetError(f'{sg_id}: no CRAM analysis, and no "reads" in assay meta')
    if meta.get('reads_type') in {'cram', 'bam'}:
        return alignment_file_rows(sg_id, meta)
    if meta.get('reads_type') in {'fastq', 'fastq_ora'}:
        return fastq_rows(sg_id, meta)
    raise SamplesheetError(f'{sg_id}: unsupported reads_type {meta.get("reads_type")!r}')


def build_rows(sgs: list[dict], access_level: str, sequencing_type: str) -> list[dict]:
    """
    Rows for every in-scope SG, each SG once, with the SG's dataset. Reports every bad SG at once, not one
    per run.
    """
    rows: list[dict] = []
    seen: set[str] = set()
    errors: list[str] = []
    for sg in sgs:
        if sg['id'] in seen:
            continue
        seen.add(sg['id'])
        # a test run must never publish to a -main bucket, nor a full run to a -test one
        project = sg['sample']['project']['name']
        if project.endswith('-test') != (access_level == 'test'):
            errors.append(f'{sg["id"]}: project {project} does not match access level {access_level}')
            continue
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type:
            logger.warning(
                f'{sg["id"]}: skipping {sg["technology"]} {sg["type"]}, this run is short-read {sequencing_type}'
            )
            continue
        try:
            rows.extend({'dataset': project.removesuffix('-test')} | row for row in rows_for(sg))
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
    access_level: str,
    sequencing_type: str,
    cohort_ids: list[str] | None,
    sg_ids: list[str] | None,
    output: str,
) -> None:
    if cohort_ids:
        sgs = fetch_cohort_sgs(sorted(set(cohort_ids)))
    else:
        sgs = fetch_sgs(sorted(set(sg_ids or [])), sequencing_type)
    rows = build_rows(sgs, access_level, sequencing_type)
    write_sheet(rows, output)
    n_sgs = len({row['sg_id'] for row in rows})
    n_datasets = len({row['dataset'] for row in rows})
    n_aligned = sum(row['kind'] == 'aligned' for row in rows)
    n_unfingerprinted = sum(row['kind'] == 'aligned' and not row['somalier'] for row in rows)
    logger.info(
        f'{sequencing_type} at access level {access_level}: {n_sgs} SGs from {n_datasets} dataset(s) '
        f'({n_aligned} already aligned, {n_unfingerprinted} of them without a somalier fingerprint), '
        f'{len(rows)} rows -> {output}',
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        '--access-level',
        required=True,
        choices=sorted(ACCESS_LEVELS),
        help='test: every SG must be in a <dataset>-test project; full or standard: none may be',
    )
    parser.add_argument('--sequencing-type', required=True, choices=sorted(SEQUENCING_TYPES))
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--cohorts', nargs='+', help='Metamist cohort IDs')
    inputs.add_argument('--sgids', nargs='+', help='Metamist sequencing group IDs, queried directly')
    parser.add_argument('--output', required=True, help='Sample sheet TSV to write')
    args = parser.parse_args()
    main(args.access_level, args.sequencing_type, args.cohorts, args.sgids, args.output)
