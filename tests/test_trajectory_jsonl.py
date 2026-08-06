"""Streaming trajectory JSONL reader/writer checks."""

import json
from pathlib import Path

import pytest
from httk.core import load

from httk.io.optimade_jsonl import TrajectoryJsonlFile, write_trajectory_jsonl


def _header(nframes: int | None = 2) -> dict[str, object]:
    return {
        "species": [{"name": "Si", "chemical_symbols": ["Si"], "concentration": [1.0]}],
        "species_at_sites": ["Si"],
        "constant_cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "nframes": nframes,
        "observable_names": ["energy"],
        "reference_frames": [0],
    }


def _frames() -> list[dict[str, object]]:
    return [
        {"fractional_site_positions": [[0.0, 0.0, 0.0]], "observables": {"energy": -1}},
        {"fractional_site_positions": [[0.5, 0.0, 0.0]], "observables": {"energy": -2}},
    ]


def test_jsonl_streams_and_registers_compression(tmp_path: Path) -> None:
    path = tmp_path / "run.traj.jsonl.gz"
    write_trajectory_jsonl(path, _header(), iter(_frames()))
    source = TrajectoryJsonlFile(path)
    assert source.path == str(path)
    assert source.header["layout"] == "dense"
    assert source.header["x-httk-trajectory"]["format"] == "httk-trajectory-jsonl"
    assert source.nframes == 2
    assert list(source.frames())[1]["index"] == 1
    assert source.frame(-1)["observables"]["energy"] == -2.0
    assert load(path, raw=True)["trajectory_jsonl"].nframes == 2


def test_jsonl_header_count_mismatch_is_an_issue(tmp_path: Path) -> None:
    path = tmp_path / "bad.jsonl"
    write_trajectory_jsonl(path, {**_header(2), "nframes": None}, _frames())
    lines = path.read_text(encoding="utf-8").splitlines()
    lines[0] = json.dumps(
        {**json.loads(lines[0]), "x-httk-trajectory": {**json.loads(lines[0])["x-httk-trajectory"], "nframes": 3}}
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    source = TrajectoryJsonlFile(path)
    assert source.nframes == 3
    assert "nframes=3" in source.issues[0]


def test_jsonl_rejects_reads_after_close(tmp_path: Path) -> None:
    path = tmp_path / "run.jsonl"
    write_trajectory_jsonl(path, _header(), _frames())
    source = TrajectoryJsonlFile(path)
    source.close()
    with pytest.raises(ValueError, match="closed"):
        _ = source.header
