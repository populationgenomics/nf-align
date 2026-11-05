"""Offline `nextflow run -stub` checks of main.nf's wiring: branching by input kind, lane/shard fan-out and
merges, CRAM QC, per-dataset published paths, registration, the run-records and pre-flight clash checks and
-resume, and the published sample sheet. The SGs span two datasets (conftest.cohort_sgs). Placeholder files stand
in for data, and the tests are skipped when nextflow isn't on PATH."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
from build_samplesheet import build_rows, write_sheet
from conftest import OTHER, PROJECT, touch

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')

MAIN_NF = Path(__file__).parents[1] / 'main.nf'
NEW = {'CPG1', 'CPG2', 'CPG3', 'CPG4', 'CPG6'}  # CPG5 and CPG7 already have a CRAM, CPG7 without a fingerprint
QC_ONLY = 'CPG7'  # already aligned, without CRAM QC: QC'd with --qc_aligned, against --aligned_ref_fasta
NEW_BY_DATASET = {PROJECT: {'CPG1', 'CPG2', 'CPG3'}, OTHER: {'CPG4', 'CPG6'}}
DATASET = {sg_id: ds for ds, sg_ids in NEW_BY_DATASET.items() for sg_id in sg_ids}


def read_tsv(path: Path, sep: str = '\t') -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


# each dataset's outputs are under <output_root>/<ds>, and the run records under the analysis dataset's
def cram_dir(out: Path, ds: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'cram'


def qc_dir(out: Path, ds: str, sg_id: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'cram_qc' / sg_id


# an SG's CRAM QC files, by their secondary-file key in its cram_qc analysis
def qc_files(sg_id: str, coverage: str = 'wgs') -> dict[str, str]:
    return {
        'verifybamid': f'{sg_id}.verify-bamid.selfSM',
        'samtools_stats': f'{sg_id}.samtools-stats',
        **{m: f'{sg_id}.{m}' for m in METRICS},
        f'picard_{coverage}_metrics': f'{sg_id}.picard-{coverage}-metrics',
        **{key: f'{sg_id}.{suffix}' for key, suffix in MOSDEPTH.items()},
        'markduplicates_metrics': f'{sg_id}.markduplicates-metrics',
    }


METRICS = [
    'alignment_summary_metrics',
    'base_distribution_by_cycle_metrics',
    'insert_size_metrics',
    'quality_by_cycle_metrics',
    'quality_yield_metrics',
]
MOSDEPTH = {
    'mosdepth_global_dist': 'mosdepth.global.dist.txt',
    'mosdepth_region_dist': 'mosdepth.region.dist.txt',
    'mosdepth_summary': 'mosdepth.summary.txt',
    'mosdepth_regions': 'regions.bed.gz',
    'mosdepth_quantized': 'quantized.bed.gz',
}
# published beside mosdepth's BEDs, but not registered
MOSDEPTH_INDEXES = ['regions.bed.gz.csi', 'quantized.bed.gz.csi']
QC_PROCESSES = ['VERIFYBAMID', 'SAMTOOLS_STATS', 'PICARD_METRICS', 'MOSDEPTH']


def published(sg_id: str, files: dict[str, str]) -> set[str]:
    """The names in an SG's QC folder: its registered QC files, and mosdepth's unregistered BED indexes."""
    return set(files.values()) | {f'{sg_id}.{index}' for index in MOSDEPTH_INDEXES}


def run_dir(out: Path, ar_guid: str = 'local-test') -> Path:
    return out / PROJECT / 'alignment_samplesheets' / ar_guid


@pytest.fixture
def run(tmp_path, cohort_sgs):
    """
    Run main.nf in stub mode over a sheet with one SG per input kind, across two datasets, with PROJECT as the
    analysis dataset; returns (result, output_root, trace rows). BUILD_SAMPLESHEET's stub hands over sheet.tsv,
    which a test may edit before its first run.
    """
    refs = tmp_path / 'refs'
    sheet = tmp_path / 'sheet.tsv'
    write_sheet(build_rows(cohort_sgs, 'standard', 'genome', str(refs / 'masked.fasta')), str(sheet))
    sheet_config = tmp_path / 'sheet.config'
    sheet_config.write_text(f"process {{ withName: 'BUILD_SAMPLESHEET' {{ ext.stub_samplesheet = '{sheet}' }} }}\n")
    for name in ['ref.fasta', 'ref.fasta.fai', *(f'ref.fasta.{ext}' for ext in ['bwt.2bit.64', 'amb', 'ann', 'pac'])]:
        touch(refs / name)
    touch(refs / 'ref.dict')
    # the reference the sample sheet gives an already-aligned CRAM whose analysis doesn't record one, with its .dict
    for name in ['masked.fasta', 'masked.fasta.fai', 'masked.dict']:
        touch(refs / name)
    touch(refs / 'sites.vcf.gz')
    # CRAM QC's references, for both sequencing types
    qc_refs = [f'{seq}_contam_{ext}' for seq in ['genome', 'exome'] for ext in ['ud', 'bed', 'mu']]
    qc_refs += ['genome_coverage_interval_list', 'exome_evaluation_interval_list']
    for name in qc_refs:
        touch(refs / name)

    def _run(
        *args: str,
        sequencing_type: str = 'genome',
        inputs: tuple = ('--cohorts', 'COH1,COH2'),
        reference: tuple = ('--ref_fasta', str(refs / 'ref.fasta')),
        ar_guid: str = 'local-test',
        omit: tuple = (),  # CRAM QC references left off the command line
    ):
        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), '-stub', '-profile', 'test',
                '-with-trace', str(trace),
                '-c', str(sheet_config),
                *inputs,
                '--output_root', str(out),
                '--dataset', PROJECT,
                '--arGuid', ar_guid,
                '--sequencing_type', sequencing_type,
                '--realignment_shards', '2',
                *reference,
                '--somalier_sites', str(refs / 'sites.vcf.gz'),
                '--aligned_ref_fasta', str(refs / 'masked.fasta'),
                '--ref_fasta_dict', str(refs / 'ref.dict'),
                *(arg for name in qc_refs if name not in omit for arg in (f'--{name}', str(refs / name))),
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


def completed(tasks: list[dict]) -> set[str]:
    """Names of the tasks that ran, other than building the sample sheet, which every new run does first."""
    return {t['name'] for t in tasks if t['status'] == 'COMPLETED' and t['name'] != 'BUILD_SAMPLESHEET (standard)'}


def task_dir(tasks: list[dict], tmp_path: Path, name: str) -> Path:
    """The work dir of the task with this name, e.g. 'REGISTER_CRAM (CPG1)'."""
    workdir = next(t for t in tasks if t['name'].startswith(name))['hash']
    return next((tmp_path / 'work').glob(f'{workdir}*'))


def test_every_input_kind_is_aligned_published_and_registered(run, tmp_path):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr

    # 2 lanes (CPG1) and 2 bazam shards (CPG3, CPG4) are merged; single units align straight to CRAM
    assert task_names(tasks, 'ALIGN_SHARD') == {'CPG1:1', 'CPG1:2', 'CPG3:1', 'CPG3:2', 'CPG4:1', 'CPG4:2'}
    assert task_names(tasks, 'MERGE_TO_CRAM') == {'CPG1', 'CPG3', 'CPG4'}
    assert task_names(tasks, 'ALIGN_TO_CRAM') == {'CPG2', 'CPG6'}
    assert task_names(tasks, 'SORT_INDEX') == {'CPG4'}  # the unindexed BAM
    # each new CRAM is fingerprinted and QC'd, then registered once, and its QC once; the aligned SGs get no
    # tasks, fingerprint or not
    assert task_names(tasks, 'SOMALIER_EXTRACT') == NEW
    for process in QC_PROCESSES:
        assert task_names(tasks, process) == NEW
    assert task_names(tasks, 'REGISTER_CRAM') == NEW
    assert task_names(tasks, 'REGISTER_CRAM_QC') == NEW

    # each SG goes to its own dataset, its somalier fingerprint beside the CRAM, and all its QC in one folder
    for ds, sg_ids in NEW_BY_DATASET.items():
        crams = {f'{s}.cram{ext}' for s in sg_ids for ext in ['', '.crai', '.somalier']}
        assert {p.name for p in cram_dir(out, ds).iterdir()} == crams
        assert {p.name for p in (out / ds / 'cram_qc').iterdir()} == sg_ids
        for sg_id in sg_ids:
            assert {p.name for p in qc_dir(out, ds, sg_id).iterdir()} == published(sg_id, qc_files(sg_id))
    assert not (out / PROJECT / 'qc').exists()

    # the run records, just the sample sheet, go to the analysis dataset only
    assert {p.name for p in run_dir(out).iterdir()} == {'samplesheet.tsv'}
    assert not (out / OTHER / 'alignment_samplesheets').exists()
    assert (run_dir(out) / 'samplesheet.tsv').read_text() == (tmp_path / 'sheet.tsv').read_text()


def test_registration_uses_the_published_path_and_the_sgs_project(run, tmp_path):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    # one analysis per SG: the CRAM, with its index and fingerprint as secondary files
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_CRAM ({sg_id})') / '.command.out').read_text()
        cram = cram_dir(out, ds) / f'{sg_id}.cram'
        assert f'--output {cram} --type cram --secondary index={cram}.crai somalier={cram}.somalier ' in command
        assert f'--sgs {sg_id} ' in command
        assert f'--project {ds} ' in command
        assert 'ar-guid=local-test' in command
        assert f'ref_fasta={tmp_path}/refs/ref.fasta' in command


def test_cram_qc_is_registered_as_one_analysis_of_the_cram(run, tmp_path):
    result, out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    # the CRAM is the primary file, and every QC file a secondary file
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_CRAM_QC ({sg_id})') / '.command.out').read_text()
        folder = qc_dir(out, ds, sg_id)
        secondary = ' '.join(f'{key}={folder / name}' for key, name in qc_files(sg_id).items())
        assert f'--output {cram_dir(out, ds) / f"{sg_id}.cram"} --type cram_qc --secondary {secondary} ' in command
        assert f'--project {ds} --' in command
        assert command.rstrip().endswith(f'--sgs {sg_id} --meta sequencing_type=genome ar-guid=local-test')


def test_cram_qc_reads_the_alignment_reference_and_the_sequencing_types_resources(run, tmp_path):
    result, _out, tasks = run(sequencing_type='exome')
    assert result.returncode == 0, result.stdout + result.stderr
    for process in QC_PROCESSES:
        staged = task_dir(tasks, tmp_path, f'CRAM_QC:{process} (CPG3)') / 'ref'
        assert {p.name for p in staged.iterdir()} == {'ref.fasta', 'ref.fasta.fai', 'ref.dict'}
    verify = task_dir(tasks, tmp_path, 'CRAM_QC:VERIFYBAMID (CPG3)')
    assert {'exome_contam_ud', 'exome_contam_bed', 'exome_contam_mu'} <= {p.name for p in verify.iterdir()}
    picard = task_dir(tasks, tmp_path, 'CRAM_QC:PICARD_METRICS (CPG3)') / 'intervals'
    assert [p.name for p in picard.iterdir()] == ['exome_evaluation_interval_list']


def test_cram_qc_can_be_turned_off(run):
    result, out, tasks = run('--cram_qc', 'false', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    for process in [*QC_PROCESSES, 'REGISTER_CRAM_QC']:
        assert not task_names(tasks, process)
    assert task_names(tasks, 'REGISTER_CRAM') == NEW
    # the duplicate metrics come from alignment, so they're still published, alone in the QC folder
    assert [p.name for p in qc_dir(out, PROJECT, 'CPG1').iterdir()] == ['CPG1.markduplicates-metrics']


def qc_tool_files(sg_id: str, coverage: str = 'wgs') -> dict[str, str]:
    """A QC-only SG's QC files: its CRAM QC tools' outputs, without the duplicate metrics alignment made."""
    return {key: name for key, name in qc_files(sg_id, coverage).items() if key != 'markduplicates_metrics'}


def existing_cram(tmp_path: Path, sg_id: str) -> Path:
    """An already-aligned SG's registered CRAM (conftest.cohort_sgs)."""
    return tmp_path / 'data' / 'existing' / f'{sg_id}.cram'


def test_aligned_sgs_without_cram_qc_are_qcd_without_realigning(run, tmp_path):
    result, out, tasks = run('--qc_aligned', 'true', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    # CPG7 is QC'd, and its QC registered, but it isn't realigned, fingerprinted or registered again. CPG5's CRAM
    # has CRAM QC already, so it gets nothing
    for process in [*QC_PROCESSES, 'REGISTER_CRAM_QC']:
        assert task_names(tasks, process) == NEW | {QC_ONLY}
    for process in ['ALIGN_TO_CRAM', 'MERGE_TO_CRAM', 'SOMALIER_EXTRACT', 'REGISTER_CRAM']:
        assert not {name.split(':')[0] for name in task_names(tasks, process)} & {QC_ONLY, 'CPG5'}

    # QC'd against the reference its CRAM was written with, which its analysis doesn't record
    for process in QC_PROCESSES:
        staged = task_dir(tasks, tmp_path, f'CRAM_QC:{process} ({QC_ONLY})') / 'ref'
        assert {p.name for p in staged.iterdir()} == {'masked.fasta', 'masked.fasta.fai', 'masked.dict'}

    # only its QC tools' files are published; its CRAM is the registered one, and there are no duplicate metrics
    assert {p.name for p in qc_dir(out, PROJECT, QC_ONLY).iterdir()} == published(QC_ONLY, qc_tool_files(QC_ONLY))
    assert not (cram_dir(out, PROJECT) / f'{QC_ONLY}.cram').exists()
    assert not (out / OTHER / 'cram_qc' / 'CPG5').exists()
    command = (task_dir(tasks, tmp_path, f'REGISTER_CRAM_QC ({QC_ONLY})') / '.command.out').read_text()
    folder = qc_dir(out, PROJECT, QC_ONLY)
    secondary = ' '.join(f'{key}={folder / name}' for key, name in qc_tool_files(QC_ONLY).items())
    assert f'--output {existing_cram(tmp_path, QC_ONLY)} --type cram_qc --secondary {secondary} --sgs' in command
    assert f'--project {PROJECT} --' in command


def test_a_run_with_nothing_to_align_only_qcs(run, tmp_path):
    sheet = tmp_path / 'sheet.tsv'
    write_sheet([r for r in read_tsv(sheet) if r['kind'] == 'aligned'], str(sheet))
    result, _out, tasks = run('--qc_aligned', 'true', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    processes = {t['name'].split(' ')[0].split(':')[-1] for t in tasks}
    assert processes == {'BUILD_SAMPLESHEET', 'PREFLIGHT', *QC_PROCESSES, 'REGISTER_CRAM_QC'}
    assert task_names(tasks, 'REGISTER_CRAM_QC') == {QC_ONLY}


@pytest.mark.parametrize('rel', ['cram_qc/CPG7', 'qc/markduplicates_metrics'])
def test_qc_only_registration_takes_the_duplicate_metrics_its_aligning_run_left(run, tmp_path, rel):
    # nf-align's own folder or cpg-flow's: not this run's file, so neither a clash nor republished
    metrics = Path(touch(tmp_path / 'results' / PROJECT / rel / f'{QC_ONLY}.markduplicates-metrics'))
    result, out, tasks = run('--qc_aligned', 'true', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    command = (task_dir(tasks, tmp_path, f'REGISTER_CRAM_QC ({QC_ONLY})') / '.command.out').read_text()
    folder = qc_dir(out, PROJECT, QC_ONLY)
    secondary = ' '.join(f'{key}={folder / name}' for key, name in qc_tool_files(QC_ONLY).items())
    assert f'--secondary {secondary} markduplicates_metrics={metrics} --sgs' in command


def test_existing_qc_files_of_a_qc_only_sg_fail_the_run(run, tmp_path):
    # e.g. an earlier attempt at its QC that failed partway
    touch(tmp_path / 'results' / PROJECT / 'cram_qc' / QC_ONLY / f'{QC_ONLY}.samtools-stats')
    result, out, tasks = run('--qc_aligned', 'true')
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert '1 output(s) already exist' in console
    assert f'{qc_dir(out, PROJECT, QC_ONLY)}/{QC_ONLY}.samtools-stats' in console
    assert not completed(tasks)


def test_failed_qc_of_a_qc_only_sg_is_not_registered(run, tmp_path):
    inject = tmp_path / 'fail.config'
    inject.write_text(
        'process {\n'
        f"    withName: 'PICARD_METRICS' {{ beforeScript = {{ meta.id == '{QC_ONLY}' ? 'exit 1' : 'true' }} }}\n"
        '}\n',
    )
    result, out, tasks = run('-c', str(inject), '--qc_aligned', 'true', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'REGISTER_CRAM_QC') == NEW
    assert not (qc_dir(out, PROJECT, QC_ONLY) / f'{QC_ONLY}.picard-wgs-metrics').exists()


def test_new_cram_is_fingerprinted_against_the_alignment_reference(run, tmp_path):
    # CPG3's input CRAM was written against old_ref.fasta, but its new CRAM against ref_fasta
    result, _out, tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr
    staged = task_dir(tasks, tmp_path, 'SOMALIER_EXTRACT (CPG3)') / 'ref'
    assert {p.name for p in staged.iterdir()} == {'ref.fasta', 'ref.fasta.fai'}


def test_test_access_level_registers_in_test_projects(run, tmp_path):
    result, out, tasks = run('--accessLevel', 'test', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    for sg_id, ds in [('CPG1', PROJECT), ('CPG4', OTHER)]:
        command = (task_dir(tasks, tmp_path, f'REGISTER_CRAM ({sg_id})') / '.command.out').read_text()
        assert f'--project {ds}-test ' in command
        # output_root decides the location, not the access level
        assert f'--output {cram_dir(out, ds) / f"{sg_id}.cram"} ' in command


def test_exomes_are_published_under_exome(run):
    result, out, _tasks = run(sequencing_type='exome')
    assert result.returncode == 0, result.stdout + result.stderr
    for ds, sg_ids in NEW_BY_DATASET.items():
        assert {p.name for p in cram_dir(out, ds, 'exome').glob('*.cram')} == {f'{s}.cram' for s in sg_ids}
        assert len(list(cram_dir(out, ds, 'exome').glob('*.somalier'))) == len(sg_ids)
        for sg_id in sg_ids:
            files = {p.name for p in qc_dir(out, ds, sg_id, 'exome').iterdir()}
            assert files == published(sg_id, qc_files(sg_id, coverage='hs'))
        assert not cram_dir(out, ds).exists()
        assert not (out / ds / 'cram_qc').exists()


def test_failed_sgs_are_dropped(run, tmp_path):
    # a stub never fails, so fail from config: one shard of CPG3, CPG2's fingerprint, CPG6's registration,
    # CPG1's samtools stats and CPG4's QC registration
    inject = tmp_path / 'fail.config'
    inject.write_text(
        'process {\n'
        "    withName: 'ALIGN_SHARD' { beforeScript = { meta.id == 'CPG3' && unit == 2 ? 'exit 1' : 'true' } }\n"
        "    withName: 'REGISTER_CRAM' { beforeScript = { meta.id == 'CPG6' ? 'exit 1' : 'true' } }\n"
        "    withName: 'SOMALIER_EXTRACT' { beforeScript = { meta.id == 'CPG2' ? 'exit 1' : 'true' } }\n"
        "    withName: 'SAMTOOLS_STATS' { beforeScript = { meta.id == 'CPG1' ? 'exit 1' : 'true' } }\n"
        "    withName: 'REGISTER_CRAM_QC' { beforeScript = { meta.id == 'CPG4' ? 'exit 1' : 'true' } }\n"
        '}\n',
    )
    result, out, tasks = run('-c', str(inject), '--register', 'true')

    # ignored failures don't fail the run
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'CPG3' not in task_names(tasks, 'MERGE_TO_CRAM')  # never merged from the surviving shard
    assert not (cram_dir(out, PROJECT) / 'CPG3.cram').exists()
    assert (cram_dir(out, OTHER) / 'CPG6.cram').exists()
    # CPG2's CRAM is published, but without a fingerprint it isn't registered, and the next run clashes on it
    assert (cram_dir(out, PROJECT) / 'CPG2.cram').exists()
    assert not (cram_dir(out, PROJECT) / 'CPG2.cram.somalier').exists()
    assert task_names(tasks, 'REGISTER_CRAM') == NEW - {'CPG2', 'CPG3'}
    # CRAM QC is all or nothing, and waits for the CRAM's own registration; a QC failure doesn't hold up the CRAM
    assert 'CPG1' in task_names(tasks, 'REGISTER_CRAM')
    assert task_names(tasks, 'REGISTER_CRAM_QC') == NEW - {'CPG1', 'CPG2', 'CPG3', 'CPG6'}
    assert (qc_dir(out, PROJECT, 'CPG1') / 'CPG1.verify-bamid.selfSM').exists()
    assert {p.name for p in run_dir(out).iterdir()} == {'samplesheet.tsv'}


def test_resume_finishes_the_sgs_that_failed(run, tmp_path):
    # CPG2's CRAM is published, but its fingerprint fails, so it isn't registered
    inject = tmp_path / 'fail.config'
    inject.write_text(
        "process { withName: 'SOMALIER_EXTRACT' { beforeScript = { meta.id == 'CPG2' ? 'exit 1' : 'true' } } }\n",
    )
    result, out, tasks = run('-c', str(inject), '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'CPG2' not in task_names(tasks, 'REGISTER_CRAM')
    assert (cram_dir(out, PROJECT) / 'CPG2.cram').exists()

    # the sample sheet is cached, so the resumed run carries on with the same SGs, though Metamist may have changed.
    # The pre-flight check passed in the first attempt and is cached, so CPG2's published CRAM isn't a clash, and
    # only its failed task and the ones after it run
    result, _out, tasks = run('-resume', '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert {'BUILD_SAMPLESHEET (standard)', 'PREFLIGHT'} <= {t['name'] for t in tasks if t['status'] == 'CACHED'}
    ran = {t['name'] for t in tasks if t['status'] == 'COMPLETED'}
    assert ran == {'SOMALIER_EXTRACT (CPG2)', 'REGISTER_CRAM (CPG2)', 'REGISTER_CRAM_QC (CPG2)'}


def test_existing_outputs_fail_the_run_before_any_alignment(run):
    _result, out, _tasks = run()
    touch(qc_dir(out, OTHER, 'CPG9') / 'CPG9.markduplicates-metrics')  # unrelated, ignored
    (cram_dir(out, OTHER) / 'CPG4.cram').unlink()  # one SG partly present still clashes on its index

    # a new arGuid, so the run-records check passes
    result, _out, tasks = run(ar_guid='second')
    assert result.returncode != 0
    # log.warn/log.error go to stderr, or to stdout when nextflow detects an AI agent (e.g. CLAUDECODE is set)
    console = result.stdout + result.stderr
    # 5 SGs x 19 outputs (2 of them unregistered BED indexes), less the deleted CRAM: a published but unregistered
    # CRAM is a clash, never registered
    assert '94 output(s) already exist' in console
    assert f'{cram_dir(out, OTHER)}/CPG4.cram.crai' in console
    assert f'{cram_dir(out, PROJECT)}/CPG1.cram' in console
    assert f'{cram_dir(out, PROJECT)}/CPG1.cram.somalier' in console
    assert f'{qc_dir(out, PROJECT, "CPG1")}/CPG1.markduplicates-metrics' in console
    assert f'{qc_dir(out, PROJECT, "CPG1")}/CPG1.verify-bamid.selfSM' in console
    assert 'CPG9' not in console
    assert not completed(tasks)


def test_resume_checks_again_after_a_failed_pre_flight_check(run, tmp_path):
    # a failed check isn't cached, so a resumed run still stops on a file it didn't write
    clash = Path(touch(cram_dir(tmp_path / 'results', PROJECT) / 'CPG1.cram'))
    for args in [(), ('-resume',)]:
        result, _out, tasks = run(*args)
        assert result.returncode != 0
        assert '1 output(s) already exist' in result.stdout + result.stderr
        assert str(clash) in result.stdout + result.stderr
        assert not completed(tasks)


def test_reusing_an_ar_guid_fails_before_any_task(run):
    result, out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    result, _out, tasks = run()
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert 'Run records for arGuid local-test already exist' in console
    assert f'{run_dir(out)}/samplesheet.tsv' in console
    assert not tasks


def test_resume_skips_the_pre_flight_check_it_passed(run):
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    # the outputs and run records exist, but the first attempt passed the pre-flight check, which is cached, and
    # -resume skips the run-records check
    result, out, tasks = run('-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert {'BUILD_SAMPLESHEET (standard)', 'PREFLIGHT'} <= {t['name'] for t in tasks}
    assert {t['status'] for t in tasks} == {'CACHED'}
    for ds, sg_ids in NEW_BY_DATASET.items():
        assert {p.name for p in cram_dir(out, ds).glob('*.cram')} == {f'{s}.cram' for s in sg_ids}


@pytest.mark.parametrize('inputs', [('--cohorts', 'COH1,COH2'), ('--sgids', 'CPG1,CPG2')])
def test_cohorts_or_sgids_are_built_into_the_published_sheet(run, inputs):
    result, out, tasks = run(inputs=inputs)
    assert result.returncode == 0, result.stdout + result.stderr
    assert task_names(tasks, 'BUILD_SAMPLESHEET') == {'standard'}
    assert (run_dir(out) / 'samplesheet.tsv').exists()


@pytest.mark.parametrize('inputs', [(), ('--cohorts', 'COH1', '--sgids', 'CPG1')])
def test_exactly_one_of_cohorts_or_sgids_is_needed(run, inputs):
    result, _out, tasks = run(inputs=inputs)
    assert result.returncode != 0
    assert 'Provide exactly one of --cohorts or --sgids' in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize(
    ('args', 'message'),
    [
        (('--sequencing_type', 'wgs'), "--sequencing_type must be genome or exome, got 'wgs'"),
        (('--accessLevel', 'prod'), "--accessLevel must be full, standard or test, got 'prod'"),
        (('--dupblaster_extra_args',), '--dupblaster_extra_args must be a string of dupblaster arguments'),
        (('--cram_qc', 'maybe'), 'Parameter `cram_qc` with type Boolean cannot be assigned to maybe'),
        (('--qc_aligned', 'maybe'), 'Parameter `qc_aligned` with type Boolean cannot be assigned to maybe'),
        (('--qc_aligned', 'true', '--cram_qc', 'false'), "--qc_aligned needs CRAM QC, so it can't be combined"),
    ],
)
def test_invalid_params_fail_early(run, args, message):
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize(('kind', 'sg_id'), [('cram', 'CPG3'), ('aligned', 'CPG5')])
def test_cram_row_without_a_reference_fails_early(run, tmp_path, kind, sg_id):
    # the sample sheet builder fills it in for an aligned CRAM whose analysis doesn't record one
    sheet = tmp_path / 'sheet.tsv'
    rows = read_tsv(sheet)
    for row in rows:
        if row['sg_id'] == sg_id:
            row['reads_ref'] = ''
    write_sheet(rows, str(sheet))
    result, _out, tasks = run()
    assert result.returncode != 0
    assert f'{sg_id}: {kind} rows need reads_ref' in result.stdout + result.stderr
    assert not completed(tasks)


def test_row_without_a_dataset_fails_early(run, tmp_path):
    sheet = tmp_path / 'sheet.tsv'
    rows = read_tsv(sheet)
    rows[0]['dataset'] = ''
    write_sheet(rows, str(sheet))
    result, _out, tasks = run()
    assert result.returncode != 0
    assert 'CPG1: a row needs dataset, the dataset the SG belongs to' in result.stdout + result.stderr
    assert not completed(tasks)


def test_per_sg_tasks_are_labelled_with_the_sgs_dataset(run, tmp_path):
    inject = tmp_path / 'labels.config'
    inject.write_text(
        'process {\n'
        "    withName: 'SORT_INDEX|ALIGN_TO_CRAM|ALIGN_SHARD|MERGE_TO_CRAM|SOMALIER_EXTRACT|VERIFYBAMID|SAMTOOLS_STATS|"
        "PICARD_METRICS|MOSDEPTH|REGISTER_CRAM|REGISTER_CRAM_QC' {\n"
        '        beforeScript = { "echo LABEL=${task.resourceLabels.dataset}" }\n'
        '    }\n'
        '}\n',
    )
    result, _out, tasks = run('-c', str(inject), '--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    labelled = [t for t in tasks if t['name'] not in {'BUILD_SAMPLESHEET (standard)', 'PREFLIGHT'}]
    # SORT_INDEX 1, ALIGN_SHARD 6, MERGE_TO_CRAM 3, ALIGN_TO_CRAM 2, then 5 each of SOMALIER_EXTRACT, the four
    # QC processes, REGISTER_CRAM and REGISTER_CRAM_QC
    assert len(labelled) == 47
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
    assert 'Parameter `cpg_container` is required' in result.stdout + result.stderr
    assert not tasks


def test_cram_qc_needs_its_containers_and_references(run, tmp_path):
    # the test profile sets placeholder containers, so unset one from config
    unset = tmp_path / 'no_picard.config'
    unset.write_text('params.picard_container = null\n')
    result, _out, tasks = run('-c', str(unset))
    assert result.returncode != 0
    assert 'CRAM QC needs --picard_container (or --cram_qc false)' in result.stdout + result.stderr
    assert not tasks

    # the references have gs:// defaults, so unset one from config and leave it off the command line
    unset.write_text('params.genome_contam_mu = null\n')
    result, _out, tasks = run('-c', str(unset), omit=('genome_contam_mu',))
    assert result.returncode != 0
    assert 'CRAM QC needs --genome_contam_mu (or --cram_qc false)' in result.stdout + result.stderr
    assert not tasks
    # an exome run doesn't need the genome references, and only CRAM QC needs any of them
    for args, seq in [((), 'exome'), (('--cram_qc', 'false'), 'genome')]:
        omit = ('genome_contam_mu',)
        result, _out, _tasks = run('-c', str(unset), *args, sequencing_type=seq, omit=omit, ar_guid=seq)
        assert result.returncode == 0, result.stdout + result.stderr


def test_missing_required_params_fail_early(run, tmp_path):
    # ref_fasta has a default, so unset it from config
    unset = tmp_path / 'no_ref.config'
    unset.write_text('params.ref_fasta = null\n')
    result, _out, tasks = run('-c', str(unset), reference=())
    assert result.returncode != 0
    assert 'Parameter `ref_fasta` is required' in result.stdout + result.stderr
    assert not tasks
