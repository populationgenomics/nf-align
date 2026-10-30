"""Offline `nextflow run -stub` checks of main.nf's wiring: branching by input kind, lane/shard fan-out and
merges, per-dataset published paths, registration, the run-records and pre-flight clash checks and -resume,
and the run records. The SGs span two datasets (conftest.cohort_sgs). Placeholder files stand in for data,
and the tests are skipped when nextflow isn't on PATH."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
from build_samplesheet import build_rows, write_sheet
from conftest import OTHER, PROJECT, touch

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')

MAIN_NF = Path(__file__).parents[1] / 'main.nf'
NEW = {'CPG1', 'CPG2', 'CPG3', 'CPG4', 'CPG6'}  # CPG5 already has a CRAM and a fingerprint
NEW_BY_DATASET = {PROJECT: {'CPG1', 'CPG2', 'CPG3'}, OTHER: {'CPG4', 'CPG6'}}
GAPFILL = 'CPG7'  # already has a CRAM, but no fingerprint
DATASET = {sg_id: ds for ds, sg_ids in NEW_BY_DATASET.items() for sg_id in sg_ids} | {GAPFILL: PROJECT}


def read_tsv(path: Path, sep: str = '\t') -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


# each dataset's outputs are under <output_root>/<ds>, and the run records under the analysis dataset's
def cram_dir(out: Path, ds: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'cram'


def metrics_dir(out: Path, ds: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'qc' / 'markduplicates_metrics'


def run_dir(out: Path, ar_guid: str = 'local-test') -> Path:
    return out / PROJECT / 'alignment_samplesheets' / ar_guid


@pytest.fixture
def run(tmp_path, cohort_sgs):
    """
    Run main.nf in stub mode over a sheet with one SG per input kind, across two datasets, with PROJECT as the
    analysis dataset; returns (result, output_root, trace rows).
    """
    sheet = tmp_path / 'sheet.tsv'
    write_sheet(build_rows(cohort_sgs, 'standard', 'genome'), str(sheet))
    refs = tmp_path / 'refs'
    for name in ['ref.fasta', 'ref.fasta.fai', *(f'ref.fasta.{ext}' for ext in ['bwt.2bit.64', 'amb', 'ann', 'pac'])]:
        touch(refs / name)
    touch(refs / 'sites.vcf.gz')
    # for the gapfilled CRAM, whose analysis doesn't record its reference
    touch(refs / 'masked.fasta')
    touch(refs / 'masked.fasta.fai')

    def _run(
        *args: str,
        sequencing_type: str = 'genome',
        inputs: tuple = ('--samplesheet', str(sheet)),
        reference: tuple = ('--ref_fasta', str(refs / 'ref.fasta')),
        profile: str = 'test',
        ar_guid: str = 'local-test',
    ):
        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), '-stub', '-profile', profile,
                '-with-trace', str(trace),
                *inputs,
                '--output_root', str(out),
                '--dataset', PROJECT,
                '--arGuid', ar_guid,
                '--sequencing_type', sequencing_type,
                '--realignment_shards', '2',
                *reference,
                '--somalier_sites', str(refs / 'sites.vcf.gz'),
                '--aligned_ref_fasta', str(refs / 'masked.fasta'),
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


def task_dir(tasks: list[dict], tmp_path: Path, name: str) -> Path:
    """The work dir of the task with this name, e.g. 'REGISTER_ANALYSIS (CPG1)'."""
    workdir = next(t for t in tasks if t['name'].startswith(name))['hash']
    return next((tmp_path / 'work').glob(f'{workdir}*'))


def test_every_input_kind_is_aligned_published_and_registered(run):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr

    # 2 lanes (CPG1) and 2 bazam shards (CPG3, CPG4) are merged; single units align straight to CRAM
    assert task_names(tasks, 'ALIGN_SHARD') == {'CPG1:1', 'CPG1:2', 'CPG3:1', 'CPG3:2', 'CPG4:1', 'CPG4:2'}
    assert task_names(tasks, 'MERGE_TO_CRAM') == {'CPG1', 'CPG3', 'CPG4'}
    assert task_names(tasks, 'ALIGN_TO_CRAM') == {'CPG2', 'CPG6'}
    assert task_names(tasks, 'SORT_INDEX') == {'CPG4'}  # the unindexed BAM
    # the gapfilled SG is only fingerprinted
    assert task_names(tasks, 'SOMALIER_EXTRACT') == NEW | {GAPFILL}
    assert task_names(tasks, 'REGISTER_ANALYSIS') == (
        {f'{s}:cram' for s in NEW} | {f'{s}:somalier' for s in NEW | {GAPFILL}}
    )

    # each SG goes to its own dataset, its somalier fingerprint beside the CRAM
    for ds, sg_ids in NEW_BY_DATASET.items():
        crams = {f'{s}.cram{ext}' for s in sg_ids for ext in ['', '.crai', '.somalier']}
        if ds == PROJECT:
            crams.add(f'{GAPFILL}.cram.somalier')
        assert {p.name for p in cram_dir(out, ds).iterdir()} == crams
        assert {p.name for p in metrics_dir(out, ds).iterdir()} == {f'{s}.markduplicates-metrics' for s in sg_ids}

    # the run records go to the analysis dataset only
    index = read_tsv(run_dir(out) / 'crams.csv', sep=',')
    assert list(index[0]) == ['sg_id', 'dataset', 'cram', 'crai']
    assert {(r['sg_id'], r['dataset']) for r in index} == {(s, DATASET[s]) for s in NEW}
    for r in index:
        assert r['cram'] == str(cram_dir(out, r['dataset']) / f'{r["sg_id"]}.cram')
        assert r['crai'] == f'{r["cram"]}.crai'
    assert (run_dir(out) / 'failed.csv').read_text() == 'sg_id,dataset,failed_at\n'
    assert not (out / OTHER / 'alignment_samplesheets').exists()
    # a sheet passed with --samplesheet isn't republished
    assert not (run_dir(out) / 'samplesheet.tsv').exists()


def test_registration_uses_the_published_path_and_the_sgs_project(run, tmp_path):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_ANALYSIS ({sg_id}:cram)') / '.command.out').read_text()
        cram = cram_dir(out, ds) / f'{sg_id}.cram'
        assert f'--output {cram} ' in command
        assert f'--type cram --secondary index={cram}.crai ' in command
        assert f'--project {ds} ' in command
        assert 'ar-guid=local-test' in command
        assert f'ref_fasta={tmp_path}/refs/ref.fasta' in command

    # a fingerprint has no index, and records the reference its CRAM was decoded with
    for sg_id, ref in [('CPG1', 'ref.fasta'), (GAPFILL, 'masked.fasta')]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_ANALYSIS ({sg_id}:somalier)') / '.command.out').read_text()
        assert f'--output {cram_dir(out, PROJECT) / f"{sg_id}.cram.somalier"} --type somalier --sgs {sg_id} ' in command
        assert f'ref_fasta={tmp_path}/refs/{ref}' in command


def test_aligned_cram_is_fingerprinted_against_its_own_reference(run, tmp_path):
    result, _out, tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr
    for sg_id, ref in [('CPG3', 'ref.fasta'), (GAPFILL, 'masked.fasta')]:
        staged = task_dir(tasks, tmp_path, f'SOMALIER_EXTRACT ({sg_id})') / 'ref'
        assert {p.name for p in staged.iterdir()} == {ref, f'{ref}.fai'}


def test_test_access_level_registers_in_test_projects(run, tmp_path):
    result, out, tasks = run('--accessLevel', 'test', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_ANALYSIS ({sg_id}:cram)') / '.command.out').read_text()
        assert f'--project {ds}-test ' in command
        # output_root decides the location, not the access level
        assert f'--output {cram_dir(out, ds) / f"{sg_id}.cram"} ' in command


def test_exomes_are_published_under_exome(run):
    result, out, _tasks = run(sequencing_type='exome')
    assert result.returncode == 0, result.stdout + result.stderr
    for ds, sg_ids in NEW_BY_DATASET.items():
        assert {p.name for p in cram_dir(out, ds, 'exome').glob('*.cram')} == {f'{s}.cram' for s in sg_ids}
        assert len(list(cram_dir(out, ds, 'exome').glob('*.somalier'))) == len(sg_ids) + (ds == PROJECT)
        assert len(list(metrics_dir(out, ds, 'exome').iterdir())) == len(sg_ids)
        assert not cram_dir(out, ds).exists()


def test_failed_sgs_are_dropped_and_reported(run, tmp_path):
    # a stub never fails, so fail from config: one shard of CPG3, CPG6's registrations, CPG2's fingerprint and
    # the gapfilled SG's fingerprint registration
    inject = tmp_path / 'fail.config'
    inject.write_text(
        'process {\n'
        "    withName: 'ALIGN_SHARD' { beforeScript = { meta.id == 'CPG3' && unit == 2 ? 'exit 1' : 'true' } }\n"
        "    withName: 'REGISTER_ANALYSIS' {\n"
        f"        beforeScript = {{ meta.id == 'CPG6' || meta.id == '{GAPFILL}' ? 'exit 1' : 'true' }}\n"
        '    }\n'
        "    withName: 'SOMALIER_EXTRACT' { beforeScript = { meta.id == 'CPG2' ? 'exit 1' : 'true' } }\n"
        '}\n',
    )
    result, out, tasks = run('-c', str(inject), '--register', 'true')

    # ignored failures don't fail the run; they're listed instead, a CRAM and its fingerprint each on their own
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'CPG3' not in task_names(tasks, 'MERGE_TO_CRAM')  # never merged from the surviving shard
    assert not (cram_dir(out, PROJECT) / 'CPG3.cram').exists()
    assert (cram_dir(out, OTHER) / 'CPG6.cram').exists()
    failed = read_tsv(run_dir(out) / 'failed.csv', sep=',')
    assert not (cram_dir(out, PROJECT) / 'CPG2.cram.somalier').exists()
    assert [(r['sg_id'], r['dataset'], r['failed_at']) for r in failed] == [
        ('CPG2', PROJECT, 'somalier'),
        ('CPG3', PROJECT, 'align'),
        ('CPG6', OTHER, 'register'),
        ('CPG6', OTHER, 'register_somalier'),
        (GAPFILL, PROJECT, 'register_somalier'),
    ]
    assert {r['sg_id'] for r in read_tsv(run_dir(out) / 'crams.csv', sep=',')} == NEW - {'CPG3'}
    assert '4 SG(s) failed and were skipped' in result.stdout + result.stderr


def test_existing_outputs_fail_the_run_before_any_alignment(run):
    _result, out, _tasks = run()
    touch(metrics_dir(out, OTHER) / 'CPG9.markduplicates-metrics')  # unrelated, ignored
    (cram_dir(out, OTHER) / 'CPG4.cram').unlink()  # one SG partly present still clashes on its index

    # a new arGuid, so the run-records check passes
    result, _out, tasks = run(ar_guid='second')
    assert result.returncode != 0
    # log.warn/log.error go to stderr, or to stdout when nextflow detects an AI agent (e.g. CLAUDECODE is set)
    console = result.stdout + result.stderr
    # the 4 SGs with a CRAM and index count as unregistered, and clash on their fingerprint alone. CPG4 clashes on
    # its index, fingerprint and metrics, and the gapfilled SG on its fingerprint
    assert '8 output(s) already exist' in console
    assert f'{cram_dir(out, OTHER)}/CPG4.cram.crai' in console
    assert f'{metrics_dir(out, OTHER)}/CPG4.markduplicates-metrics' in console
    assert f'{cram_dir(out, PROJECT)}/CPG1.cram.somalier' in console
    assert f'{cram_dir(out, PROJECT)}/{GAPFILL}.cram.somalier' in console
    assert 'CPG1.markduplicates-metrics' not in console
    assert 'CPG9' not in console
    assert not [t for t in tasks if t['status'] == 'COMPLETED']


def test_unregistered_cram_is_fingerprinted_and_registered_not_realigned(run, tmp_path):
    # an earlier run published CPG1's CRAM and index (CPG1 has 2 FASTQ lanes), and CPG4's (a BAM), but registered
    # neither. A stale metrics file doesn't matter
    earlier = tmp_path / 'results'
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        touch(cram_dir(earlier, ds) / f'{sg_id}.cram')
        touch(cram_dir(earlier, ds) / f'{sg_id}.cram.crai')
    touch(metrics_dir(earlier, PROJECT) / 'CPG1.markduplicates-metrics')

    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'ALIGN_SHARD') == {'CPG3:1', 'CPG3:2'}
    assert task_names(tasks, 'SORT_INDEX') == set()
    assert task_names(tasks, 'SOMALIER_EXTRACT') == NEW | {GAPFILL}
    # the found CRAMs are registered in tasks of their own, as are their fingerprints
    assert task_names(tasks, 'REGISTER_ANALYSIS') == (
        {f'{s}:cram' for s in NEW} | {f'{s}:somalier' for s in NEW | {GAPFILL}}
    )
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        assert (cram_dir(out, ds) / f'{sg_id}.cram.somalier').exists()
        # fingerprinted against the alignment reference, which an earlier nf-align run wrote the CRAM with
        staged = task_dir(tasks, tmp_path, f'SOMALIER_EXTRACT ({sg_id})')
        assert {p.name for p in (staged / 'ref').iterdir()} == {'ref.fasta', 'ref.fasta.fai'}
        assert (staged / f'{sg_id}.cram').resolve() == cram_dir(out, ds) / f'{sg_id}.cram'
        command = (task_dir(tasks, tmp_path, f'REGISTER_ANALYSIS ({sg_id}:cram)') / '.command.out').read_text()
        cram = cram_dir(out, ds) / f'{sg_id}.cram'
        assert f'--output {cram} --type cram --secondary index={cram}.crai ' in command
        assert f'ref_fasta={tmp_path}/refs/ref.fasta' in command
    # the run records list only the CRAMs this run aligned
    assert {r['sg_id'] for r in read_tsv(run_dir(out) / 'crams.csv', sep=',')} == NEW - {'CPG1', 'CPG4'}
    assert (run_dir(out) / 'failed.csv').read_text() == 'sg_id,dataset,failed_at\n'

    # the fingerprints this session claimed don't stop a resume, and nothing reruns
    result, _out, tasks = run('--register', 'true', '-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert {t['status'] for t in tasks} == {'CACHED'}


def test_reusing_an_ar_guid_fails_before_any_task(run):
    result, out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    result, _out, tasks = run()
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert 'Run records for arGuid local-test already exist' in console
    assert f'{run_dir(out)}/crams.csv' in console
    assert not tasks


def test_resume_republishes_outputs_this_session_claimed(run):
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    # the outputs and run records exist, but this session's earlier attempt claimed the outputs at the
    # pre-flight check, and -resume skips the run-records check
    result, out, tasks = run('-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert any(t['name'] == 'WRITE_RUN_RECORDS' for t in tasks)
    assert {t['status'] for t in tasks} == {'CACHED'}
    for ds, sg_ids in NEW_BY_DATASET.items():
        assert {p.name for p in cram_dir(out, ds).glob('*.cram')} == {f'{s}.cram' for s in sg_ids}


def test_cohorts_mode_publishes_the_built_sheet(run):
    result, out, tasks = run(inputs=('--cohorts', 'COH1,COH2'))
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'BUILD_SAMPLESHEET') == {'standard'}
    assert (run_dir(out) / 'samplesheet.tsv').exists()


def test_sgids_mode_publishes_the_built_sheet(run):
    result, out, tasks = run(inputs=('--sgids', 'CPG1,CPG2'))
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'BUILD_SAMPLESHEET') == {'standard'}
    assert (run_dir(out) / 'samplesheet.tsv').exists()


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


def test_cram_row_without_a_reference_fails_early(run, tmp_path):
    sheet = tmp_path / 'sheet.tsv'
    rows = read_tsv(sheet)
    for row in rows:
        if row['kind'] == 'cram':
            row['reads_ref'] = ''
    write_sheet(rows, str(sheet))
    result, _out, tasks = run()
    assert result.returncode != 0
    assert 'CPG3: a cram row needs reads_ref' in result.stdout + result.stderr
    assert not [t for t in tasks if t['status'] == 'COMPLETED']


def test_row_without_a_dataset_fails_early(run, tmp_path):
    sheet = tmp_path / 'sheet.tsv'
    rows = read_tsv(sheet)
    rows[0]['dataset'] = ''
    write_sheet(rows, str(sheet))
    result, _out, tasks = run()
    assert result.returncode != 0
    assert 'CPG1: a row needs dataset, the dataset the SG belongs to' in result.stdout + result.stderr
    assert not [t for t in tasks if t['status'] == 'COMPLETED']


def test_per_sg_tasks_are_labelled_with_the_sgs_dataset(run, tmp_path):
    inject = tmp_path / 'labels.config'
    inject.write_text(
        'process {\n'
        "    withName: 'SORT_INDEX|ALIGN_TO_CRAM|ALIGN_SHARD|MERGE_TO_CRAM|SOMALIER_EXTRACT|REGISTER_ANALYSIS' {\n"
        '        beforeScript = { "echo LABEL=${task.resourceLabels.dataset}" }\n'
        '    }\n'
        '}\n',
    )
    result, _out, tasks = run('-c', str(inject), '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    labelled = [t for t in tasks if not t['name'].startswith('WRITE_RUN_RECORDS')]
    # SORT_INDEX 1, ALIGN_SHARD 6, MERGE_TO_CRAM 3, ALIGN_TO_CRAM 2, SOMALIER_EXTRACT 6, REGISTER_ANALYSIS 5 + 6
    assert len(labelled) == 29
    for t in labelled:
        sg_id = t['name'].split('(')[1].rstrip(')').split(':')[0]
        # beforeScript's output goes to the task log, not .command.out
        assert f'LABEL={DATASET[sg_id]}' in (task_dir(tasks, tmp_path, t['name']) / '.command.log').read_text()


def test_workflow_needs_every_container(run, tmp_path):
    # the test profile sets placeholder containers, and an empty value can't be passed on the command line
    unset = tmp_path / 'no_cpg_container.config'
    unset.write_text('params.cpg_container = null\n')
    result, _out, tasks = run('-c', str(unset), '--align_container', 'image:1')
    assert result.returncode != 0
    message = 'The workflow needs --align_container, --cpg_container and --somalier_container'
    assert message in result.stdout + result.stderr
    assert not tasks


def test_missing_required_params_fail_early(run, tmp_path):
    # ref_fasta has a default, so unset it from config
    unset = tmp_path / 'no_ref.config'
    unset.write_text('params.ref_fasta = null\n')
    result, _out, tasks = run('-c', str(unset), reference=())
    assert result.returncode != 0
    assert 'Missing required params: --ref_fasta' in result.stdout + result.stderr
    assert not tasks
