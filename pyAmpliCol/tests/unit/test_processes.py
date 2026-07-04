from __future__ import annotations

import subprocess
from pathlib import Path

from pyamplicol.processes import ProcessEnumerator, ProcessOptions


def test_structured_process_enumeration_for_z_plus_two_gluons() -> None:
    enumeration = ProcessEnumerator().enumerate("d d~ > z g g")

    assert enumeration.request.initial_state == ("d~", "d")
    assert enumeration.request.rest == ("z", "g", "g")
    assert len(enumeration.unique_processes) == 1
    assert len(enumeration.groups) == 1
    assert enumeration.n_records == 1

    record = enumeration.groups[0].records[0]
    assert record.process == ("d~", "d", "g", "g", "z")
    assert record.identical_factor == 2.0


def test_process_export_matches_legacy_process_list_for_small_qcd_case(
    tmp_path: Path,
) -> None:
    repo_root = Path(__file__).resolve().parents[3]
    subprocess.run(
        [
            "python3",
            str(repo_root / "process_list.py"),
            "--serial",
            "d d~ > z g g",
        ],
        cwd=tmp_path,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    legacy_text = (tmp_path / "processes.txt").read_text()

    enumerator = ProcessEnumerator(ProcessOptions(serial=True))
    enumeration = enumerator.enumerate("d d~ > z g g")
    export = tmp_path / "pyamplicol-processes.txt"
    enumerator.write_legacy_file(enumeration, export)

    assert export.read_text() == legacy_text
