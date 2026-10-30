"""Metamist-shaped fixtures shared by the sample sheet and stub-run tests. Metamist is never contacted."""

from pathlib import Path

import pytest

PROJECT = 'ds'
OTHER = 'other'  # a second dataset, so runs span two


def touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return str(path)


def sg(
    sg_id: str,
    assays: list[dict] | None = None,
    analyses: list[dict] | None = None,
    project: str = PROJECT,
    **fields: object,
) -> dict:
    """A sequencing group, shaped like the SgFields fragment of the sample sheet queries."""
    return {
        'id': sg_id,
        'type': 'genome',
        'technology': 'short-read',
        'sample': {'project': {'name': project}},
        'assays': [{'meta': meta} for meta in assays or []],
        'analyses': analyses or [],
    } | fields


def cohort(cohort_id: str, sgs: list[dict], status: str = 'active') -> dict:
    return {'id': cohort_id, 'status': status, 'sequencingGroups': sgs}


def fastq_meta(r1: str, r2: str, **extra: object) -> dict:
    return {'reads_type': 'fastq', 'reads': [{'location': r1}, {'location': r2}]} | extra


@pytest.fixture
def cohort_sgs(tmp_path: Path) -> list[dict]:
    """
    One SG per input kind, as local placeholder files so a stub run can stage them. CPG1-CPG3 are in PROJECT,
    CPG4-CPG6 in OTHER.
    """
    data = tmp_path / 'data'
    cram = touch(data / 'existing' / 'CPG5.cram')
    touch(data / 'existing' / 'CPG5.cram.crai')
    touch(data / 'old_ref.fasta.fai')
    return [
        sg('CPG1', [fastq_meta(touch(data / f'L{n}_R1.fq.gz'), touch(data / f'L{n}_R2.fq.gz')) for n in (1, 2)]),
        sg(
            'CPG2',
            [
                fastq_meta(
                    touch(data / 'R1.fq.ora'),
                    touch(data / 'R2.fq.ora'),
                    reads_type='fastq_ora',
                    ora_reference={'location': touch(data / 'ora.tar')},
                ),
            ],
        ),
        sg(
            'CPG3',
            [
                {
                    'reads_type': 'cram',
                    'reads': {
                        'location': touch(data / 'CPG3.cram'),
                        'secondaryFiles': [{'location': touch(data / 'CPG3.cram.crai')}],
                    },
                    'reference_assembly': {'location': touch(data / 'old_ref.fasta')},
                },
            ],
        ),
        sg('CPG4', [{'reads_type': 'bam', 'reads': {'location': touch(data / 'CPG4.bam')}}], project=OTHER),
        sg('CPG5', analyses=[{'id': 1, 'output': cram}], project=OTHER),
        sg('CPG6', [fastq_meta(touch(data / 'S6_R1.fq.gz'), touch(data / 'S6_R2.fq.gz'))], project=OTHER),
    ]
