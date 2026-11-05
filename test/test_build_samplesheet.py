"""bin/build_samplesheet.py: Metamist responses -> sample sheet rows. The query is validated against the
schema bundled with the metamist client; responses are hand-built dicts in the shape that query returns."""

import csv

import build_samplesheet
import pytest
from build_samplesheet import SamplesheetError, build_rows, fetch_cohort_sgs, fetch_sgs, rows_for, write_sheet
from conftest import PROJECT, cohort, fastq_meta, sg
from metamist.graphql import gql, validate


@pytest.mark.parametrize('query', [build_samplesheet.COHORTS_QUERY, build_samplesheet.SGS_QUERY])
def test_queries_match_metamist_schema(query):
    assert validate(gql(query), use_local_schema=True)


def test_rows_classify_each_input_kind(cohort_sgs):
    rows = build_rows(cohort_sgs, PROJECT, 'genome')
    assert [(r['sg_id'], r['kind'], r['n_units'], r['unit']) for r in rows] == [
        ('CPG1', 'fastq', 2, 1),
        ('CPG1', 'fastq', 2, 2),
        ('CPG2', 'fastq_ora', 1, 1),
        ('CPG3', 'cram', 1, 1),
        ('CPG4', 'bam', 1, 1),
        ('CPG5', 'aligned', 1, 1),
        ('CPG6', 'fastq', 1, 1),
    ]
    assert rows[0]['r1'].endswith('L1_R1.fq.gz')
    assert rows[2]['reads_ref'].endswith('ora.tar')
    assert rows[3]['reads_ref'].endswith('old_ref.fasta')
    assert rows[3]['reads_index'].endswith('CPG3.cram.crai')
    assert rows[4]['reads_index'] == ''  # unindexed BAM goes through SORT_INDEX
    assert rows[5]['reads_index'].endswith('CPG5.cram.crai')


def test_read_locations_are_used_verbatim():
    meta = fastq_meta('gs://cpg-ds-main-upload/R1.fq.gz', 'gs://cpg-ds-main-upload/R2.fq.gz')
    assert rows_for(sg('CPG1', [meta]))[0]['r1'] == 'gs://cpg-ds-main-upload/R1.fq.gz'


def test_latest_cram_analysis_wins():
    analyses = [{'id': 9, 'output': 'gs://b/new.cram'}, {'id': 3, 'output': 'gs://b/old.cram'}]
    assert rows_for(sg('CPG1', analyses=analyses))[0]['reads'] == 'gs://b/new.cram'


def test_cram_analysis_without_output_is_ignored():
    meta = fastq_meta('gs://b/R1.fq.gz', 'gs://b/R2.fq.gz')
    assert rows_for(sg('CPG1', [meta], analyses=[{'id': 1, 'output': None}]))[0]['kind'] == 'fastq'


def test_sgs_are_deduplicated_and_filtered_to_the_sequencing_type(caplog):
    meta = fastq_meta('gs://b/R1.fq.gz', 'gs://b/R2.fq.gz')
    # CPG1 is in two cohorts
    sgs = [
        sg('CPG1', [meta]),
        sg('CPG2', [meta], type='exome'),
        sg('CPG1', [meta]),
        sg('CPG3', [meta], technology='long-read'),
    ]
    assert [r['sg_id'] for r in build_rows(sgs, PROJECT, 'genome')] == ['CPG1']
    assert [r['sg_id'] for r in build_rows(sgs, PROJECT, 'exome')] == ['CPG2']
    assert 'CPG3: skipping long-read genome' in caplog.text


def test_sgs_from_another_project_fail():
    meta = fastq_meta('gs://b/R1.fq.gz', 'gs://b/R2.fq.gz')
    other = sg('CPG2', [meta], sample={'project': {'name': 'other'}})
    with pytest.raises(SamplesheetError, match='CPG2: belongs to project other, not ds'):
        build_rows([sg('CPG1', [meta]), other], PROJECT, 'genome')


def test_no_in_scope_sgs_fails():
    meta = fastq_meta('gs://b/R1.fq.gz', 'gs://b/R2.fq.gz')
    with pytest.raises(SamplesheetError, match='No short-read exome'):
        build_rows([sg('CPG1', [meta])], PROJECT, 'exome')


def test_bad_sgs_are_reported_together():
    sgs = [
        sg('CPG1'),
        sg('CPG2', [{'reads_type': 'fastq_ora', 'reads': [{'location': 'a'}, {'location': 'b'}]}]),
        sg('CPG3', [{'reads_type': 'cram', 'reads': {'location': 'gs://b/x.bam'}}]),
        sg('CPG4', [{'reads_type': 'fastq', 'reads': [{'location': 'a'}]}]),
    ]
    with pytest.raises(SamplesheetError) as error:
        build_rows(sgs, PROJECT, 'genome')
    assert [line.split(':')[0] for line in str(error.value).splitlines()] == ['CPG1', 'CPG2', 'CPG3', 'CPG4']


def test_cram_without_a_reference_fails():
    # decoding with ref_fasta instead fails, deep in bazam, when the CRAM was written against another reference
    meta = {'reads_type': 'cram', 'reads': {'location': 'gs://b/CPG1.cram'}}
    with pytest.raises(SamplesheetError, match=r'CPG1: "reference_assembly\.location" is required for cram assays'):
        rows_for(sg('CPG1', [meta]))


def test_fetch_cohorts_is_scoped_to_the_project(monkeypatch):
    calls = []

    def fake_query(_query, variables):
        calls.append(variables)
        return {
            'project': {
                'cohorts': [
                    cohort('COH1', [sg('CPG1')]),
                    cohort('COH2', [sg('CPG2')], status='archived'),
                    cohort('COH3', []),
                ],
            },
        }

    monkeypatch.setattr(build_samplesheet, 'query', fake_query)
    # COH4 may exist in another project, but isn't visible from this one
    with pytest.raises(SamplesheetError, match=r"missing: \['COH4'\], inactive: \['COH2'\], empty: \['COH3'\]"):
        fetch_cohort_sgs('ds-test', ['COH1', 'COH2', 'COH3', 'COH4'])
    assert calls == [{'project': 'ds-test', 'cohorts': ['COH1', 'COH2', 'COH3', 'COH4']}]


def test_fetch_cohort_sgs_flattens_cohorts(monkeypatch):
    response = {'project': {'cohorts': [cohort('COH1', [sg('CPG1'), sg('CPG2')]), cohort('COH2', [sg('CPG1')])]}}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, **_variables: response)
    assert [s['id'] for s in fetch_cohort_sgs('ds', ['COH1', 'COH2'])] == ['CPG1', 'CPG2', 'CPG1']


def test_fetch_sgs_is_scoped_to_the_project(monkeypatch):
    calls = []

    def fake_query(_query, variables):
        calls.append(variables)
        return {'project': {'sequencingGroups': [sg('CPG1'), sg('CPG2')]}}

    monkeypatch.setattr(build_samplesheet, 'query', fake_query)
    assert [s['id'] for s in fetch_sgs('ds-test', ['CPG1', 'CPG2'], 'genome')] == ['CPG1', 'CPG2']
    # CPG3 may exist in another project, but isn't visible from this one
    with pytest.raises(SamplesheetError, match=r"missing: \['CPG3'\]"):
        fetch_sgs('ds-test', ['CPG1', 'CPG2', 'CPG3'], 'genome')
    assert calls[-1] == {'project': 'ds-test', 'sgids': ['CPG1', 'CPG2', 'CPG3']}


def test_fetch_sgs_fails_on_sgs_of_the_wrong_type(monkeypatch):
    # named SGs fail rather than being skipped, unlike SGs found through a cohort
    sgs = [sg('CPG1'), sg('CPG2', type='exome'), sg('CPG3', technology='long-read')]
    response = {'project': {'sequencingGroups': sgs}}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, **_variables: response)
    with pytest.raises(
        SamplesheetError,
        match=r"not short-read genome: \['CPG2 \(short-read exome\)', 'CPG3 \(long-read genome\)'\]",
    ):
        fetch_sgs('ds', ['CPG1', 'CPG2', 'CPG3'], 'genome')


def test_write_sheet_has_every_column(tmp_path, cohort_sgs):
    path = tmp_path / 'sheet.tsv'
    write_sheet(build_rows(cohort_sgs, PROJECT, 'genome'), str(path))
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert list(rows[0]) == build_samplesheet.COLUMNS
    assert rows[4]['reads_index'] == ''
