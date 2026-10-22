"""Offline `nextflow run -stub` checks of main.nf's wiring: branching by input kind, lane/shard fan-out and
merges, published paths, registration, the pre-flight clash check and -resume, and the failed-SG report.
Placeholder files stand in for data, and the tests are skipped when nextflow isn't on PATH."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
from build_samplesheet import build_rows, write_sheet
from conftest import PROJECT, touch

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')

MAIN_NF = Path(__file__).parents[1] / 'main.nf'
NEW = {'CPG1', 'CPG2', 'CPG3', 'CPG4', 'CPG6'}  # CPG5 already has a CRAM
RUN_DIR = Path('alignment_samplesheets/local-test')


def read_tsv(path: Path, sep: str = '\t') -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


@pytest.fixture
def run(tmp_path, cohort_sgs):
    """Run main.nf in stub mode over a sheet with one SG per input kind; returns (result, outputDir, trace rows)."""
    sheet = tmp_path / 'sheet.tsv'
    write_sheet(build_rows(cohort_sgs, PROJECT, 'genome'), str(sheet))
    refs = tmp_path / 'refs'
    for name in ['ref.fasta', 'ref.fasta.fai', *(f'ref.fasta.{ext}' for ext in ['bwt.2bit.64', 'amb', 'ann', 'pac'])]:
        touch(refs / name)

    def _run(
        *args: str,
        sequencing_type: str = 'genome',
        inputs: tuple = ('--samplesheet', str(sheet)),
        index: tuple = ('--bwa_index_prefix', str(refs / 'ref.fasta')),
        profile: str = 'test',
    ):
        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), '-stub', '-profile', profile,
                '-with-trace', str(trace),
                '-output-dir', str(out),
                *inputs,
                '--dataset', PROJECT,
                '--sequencing_type', sequencing_type,
                '--realignment_shards', '2',
                '--ref_fasta', str(refs / 'ref.fasta'),
                *index,
                *args,
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )  # fmt: skip
        tasks = read_tsv(trace) if trace.exists() else []
        return result, out, tasks

    return _run


def task_names(tasks: list[dict], process: str) -> set[str]:
    """Tags of the trace's tasks for one process, e.g. {'CPG1', 'CPG3:2'}."""
    return {t['name'].split('(')[1].rstrip(')') for t in tasks if t['name'].split(' ')[0].endswith(process)}


def test_every_input_kind_is_aligned_published_and_registered(run):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr

    # 2 lanes (CPG1) and 2 bazam shards (CPG3, CPG4) are merged; single units align straight to CRAM
    assert task_names(tasks, 'ALIGN_SHARD') == {'CPG1:1', 'CPG1:2', 'CPG3:1', 'CPG3:2', 'CPG4:1', 'CPG4:2'}
    assert task_names(tasks, 'MERGE_TO_CRAM') == {'CPG1', 'CPG3', 'CPG4'}
    assert task_names(tasks, 'ALIGN_TO_CRAM') == {'CPG2', 'CPG6'}
    assert task_names(tasks, 'SORT_INDEX') == {'CPG4'}  # the unindexed BAM
    assert task_names(tasks, 'REGISTER_ANALYSIS') == NEW

    assert {p.name for p in (out / 'cram').iterdir()} == {f'{s}.cram{ext}' for s in NEW for ext in ['', '.crai']}
    assert {p.name for p in (out / 'qc' / 'markduplicates_metrics').iterdir()} == {
        f'{s}.markduplicates-metrics' for s in NEW
    }
    index = read_tsv(out / RUN_DIR / 'crams.csv', sep=',')
    assert {r['id'] for r in index} == NEW
    assert all(r['cram'] == str(out / 'cram' / f'{r["id"]}.cram') for r in index)
    assert not read_tsv(out / RUN_DIR / 'failed.csv', sep=',')
    # a sheet passed with --samplesheet isn't republished
    assert not (out / RUN_DIR / 'samplesheet.tsv').exists()


def test_registration_uses_the_published_path(run, tmp_path):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    workdir = next(t for t in tasks if t['name'].startswith('REGISTER_ANALYSIS (CPG1)'))['hash']
    command = next((tmp_path / 'work').glob(f'{workdir}*/.command.out')).read_text()
    assert f'--output {out}/cram/CPG1.cram ' in command
    assert f'--secondary index={out}/cram/CPG1.cram.crai ' in command
    assert '--project ds ' in command
    assert 'ar-guid=local-test' in command


def test_exomes_are_published_under_exome(run):
    result, out, _tasks = run(sequencing_type='exome')
    assert result.returncode == 0, result.stdout + result.stderr
    assert {p.name for p in (out / 'exome' / 'cram').glob('*.cram')} == {f'{s}.cram' for s in NEW}
    assert len(list((out / 'exome' / 'qc' / 'markduplicates_metrics').iterdir())) == len(NEW)
    assert not (out / 'cram').exists()


def test_failed_sgs_are_dropped_and_reported(run, tmp_path):
    # a stub never fails, so fail one shard of CPG3 and CPG6's registration from config
    inject = tmp_path / 'fail.config'
    inject.write_text(
        'process {\n'
        "    withName: 'ALIGN_SHARD' { beforeScript = { meta.id == 'CPG3' && unit == 2 ? 'exit 1' : 'true' } }\n"
        "    withName: 'REGISTER_ANALYSIS' { beforeScript = { meta.id == 'CPG6' ? 'exit 1' : 'true' } }\n"
        '}\n',
    )
    result, out, tasks = run('-c', str(inject), '--register', 'true')

    # ignored failures don't fail the run; they're listed instead
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'CPG3' not in task_names(tasks, 'MERGE_TO_CRAM')  # never merged from the surviving shard
    assert not (out / 'cram' / 'CPG3.cram').exists()
    assert (out / 'cram' / 'CPG6.cram').exists()
    failed = read_tsv(out / RUN_DIR / 'failed.csv', sep=',')
    assert sorted((r['id'], r['failed_at']) for r in failed) == [('CPG3', 'align'), ('CPG6', 'register')]
    assert '2 SG(s) failed and were skipped' in result.stdout


def test_existing_outputs_fail_the_run_before_any_alignment(run):
    _result, out, _tasks = run()
    touch(out / 'qc' / 'markduplicates_metrics' / 'CPG9.markduplicates-metrics')  # unrelated, ignored
    (out / 'cram' / 'CPG2.cram').unlink()  # one SG partly present still clashes on its index

    result, _out, tasks = run()
    assert result.returncode != 0
    assert '14 output(s) already exist' in result.stdout
    assert f'{out}/cram/CPG2.cram.crai' in result.stdout
    assert 'CPG9' not in result.stdout
    assert not [t for t in tasks if t['status'] == 'COMPLETED']


def test_resume_republishes_outputs_this_session_claimed(run):
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    # the outputs exist, but this session's earlier attempt claimed them at the pre-flight check
    result, out, tasks = run('-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert tasks
    assert {t['status'] for t in tasks} == {'CACHED'}
    assert {p.name for p in (out / 'cram').glob('*.cram')} == {f'{s}.cram' for s in NEW}


def test_cohorts_mode_publishes_the_built_sheet(run):
    result, out, tasks = run(inputs=('--cohorts', 'COH1,COH2'))
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'BUILD_SAMPLESHEET') == {'ds'}
    assert (out / RUN_DIR / 'samplesheet.tsv').exists()


def test_sgids_mode_publishes_the_built_sheet(run):
    result, out, tasks = run(inputs=('--sgids', 'CPG1,CPG2'))
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'BUILD_SAMPLESHEET') == {'ds'}
    assert (out / RUN_DIR / 'samplesheet.tsv').exists()


@pytest.mark.parametrize(
    ('args', 'message'),
    [
        (('--sequencing_type', 'wgs'), "--sequencing_type must be genome or exome, got 'wgs'"),
        (('--accessLevel', 'prod'), "--accessLevel must be full, standard or test, got 'prod'"),
        (('--cohorts', 'COH1'), 'Provide exactly one of --cohorts, --sgids or --samplesheet'),
        (('--sgids', 'CPG1'), 'Provide exactly one of --cohorts, --sgids or --samplesheet'),
        (('--dupblaster_extra_args',), '--dupblaster_extra_args must be a string of dupblaster arguments'),
    ],
)
def test_invalid_params_fail_early(run, args, message):
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


def test_workflow_needs_both_containers(run, tmp_path):
    # the test profile sets placeholder containers, and an empty value can't be passed on the command line
    unset = tmp_path / 'no_cpg_container.config'
    unset.write_text('params.cpg_container = null\n')
    result, _out, tasks = run('-c', str(unset), '--align_container', 'image:1')
    assert result.returncode != 0
    assert 'The workflow needs --align_container and --cpg_container' in result.stdout + result.stderr
    assert not tasks


def test_missing_required_params_fail_early(run):
    result, _out, tasks = run(index=())
    assert result.returncode != 0
    assert 'Missing required params: --bwa_index_prefix' in result.stdout + result.stderr
    assert not tasks
