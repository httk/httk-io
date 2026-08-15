"""Tests for the lazy VASP OUTCAR reader."""

import bz2
import io
from pathlib import Path

import httk.core
import pytest

from httk.io.vasp import ElasticModuliBlock, OutcarFile, read_outcar

OUTCAR = """ vasp.5.2.12 01Jan20 (build synthetic) complex
   LEXCH  = PE
   TITEL  = PAW_PBE Na_pv 05Jan2001
   TITEL  = PAW_PBE Na_pv 05Jan2001
   ENCUT  = 300.0 eV
   k-points           NKPTS =     2   k-points in BZ
   POTIM  = 0.5000    time-step for ionic-motion
   NSW    =    10    number of steps for IOM
   IBRION =      2    ionic relax
   ISIF   =      3    stress and relaxation
   NELM   =     60;   NELMIN=  2
   EDIFF  = 0.1E-05   stopping-criterion for ELM
   EDIFFG = 0.1E-04   stopping-criterion for IOM
   ISMEAR =     0;   SIGMA  =   0.10  broadening
   ISPIN  =      2    spin polarized calculation?
   GGA     =    --    GGA type
 ----------------------------------------- Iteration    1(   1)  ---------------------------------------
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -27.09328752 eV
   energy  without entropy=      -27.09328752  energy(sigma->0) =      -27.09328752
  General timing and accounting informations for this job:
"""

FRAME_HEADER = """ vasp.5.2.12 synthetic complex
   LEXCH  = PE
   GGA     =    --
"""


def frame_text(number: int, energy: str, stress: tuple[str, ...], *, temperature: str | None = None) -> str:
    temperature_line = f"  kin. lattice  EKIN_LAT= 0.000000 (temperature {temperature} K)\n" if temperature else ""
    return f""" ----------------------------------------- Iteration {number:4d}(   1)  ---------------------------------------
  in kB       {' '.join(stress)}
       direct lattice vectors                 reciprocal lattice vectors
     {number}.000000000  0.000000000  0.000000000     0.1 0.0 0.0
     0.000000000  {number}.000000000  0.000000000     0.0 0.1 0.0
     0.000000000  0.000000000  {number}.000000000     0.0 0.0 0.1
  POSITION                                       TOTAL-FORCE (eV/Angst)
  -----------------------------------------------------------------------------------
     0.{number}0000  0.000000  0.000000       -0.{number}00000  0.000000  0.000000
     0.000000  0.{number}0000  0.000000        0.000000 -0.{number}00000  0.000000
    total drift:                              0.000000 0.000000 0.000000
  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  free  energy   TOTEN  =       {energy} eV
  energy  without entropy=      {energy}  energy(sigma->0) =      {energy}
{temperature_line}"""


FRAME_OUTCAR = (
    FRAME_HEADER
    + frame_text(1, "-1.100", ("1.0", "2.0", "3.0", "4.0", "5.0", "6.0"))
    + frame_text(2, "-2.200", ("2.0", "3.0", "4.0", "5.0", "6.0", "7.0"))
    + frame_text(3, "-3.300", ("3.0", "4.0", "5.0", "6.0", "7.0", "8.0"))
    + " General timing and accounting informations for this job:\n"
)

NIONS_FRAME_HEADER = FRAME_HEADER + "   number of ions     NIONS =      2\n   ions per type =               1   1\n"


def test_outcar_prologue_and_final_results(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(OUTCAR, encoding="utf-8")
    outcar = OutcarFile(path)
    assert outcar.path == str(path)
    assert outcar.version_string == "vasp.5.2.12"
    assert outcar.version_numbers == (5, 2, 12)
    assert outcar.parameters["ENCUT"] == "300.0"
    assert outcar.parameters["NKPTS"] == "2"
    assert outcar.parameters["SIGMA"] == "0.10"
    assert outcar.xc == "Perdew-Burke-Ernzerhof (PBE)"
    assert outcar.potcar_titles == ("PAW_PBE Na_pv 05Jan2001",)
    assert outcar.final_energies.free_energy == "-27.09328752"
    assert outcar.final_energies.energy_without_entropy == "-27.09328752"
    assert outcar.final_energies.energy_sigma0 == "-27.09328752"
    assert outcar.final_energies.final is True
    assert outcar.completed is True
    assert outcar.completion_evidence == ("  General timing and accounting informations for this job:",)
    assert outcar.issues == ()


def test_outcar_version_suffix_and_truncation(tmp_path: Path) -> None:
    path = tmp_path / "truncated.outcar"
    path.write_text(
        OUTCAR.replace("vasp.5.2.12", "vasp.5.4.4.18Apr17-6-g9f103f2a35").split("  General timing")[0], encoding="utf-8"
    )
    outcar = OutcarFile(path)
    assert outcar.version_string == "vasp.5.4.4.18Apr17-6-g9f103f2a35"
    assert outcar.version_numbers == (5, 4, 4)
    assert outcar.final_energies.final is True
    assert outcar.completed is False
    assert any("missing completion footer" in issue for issue in outcar.issues)


def test_outcar_unknown_xc_and_unparseable_parameter(tmp_path: Path) -> None:
    path = tmp_path / "unknown.outcar"
    path.write_text(
        OUTCAR.replace("GGA     =    --", "GGA     =    UNKNOWN").replace("ENCUT  = 300.0", "ENCUT  = ???"),
        encoding="utf-8",
    )
    outcar = OutcarFile(path)
    assert outcar.xc is None
    assert outcar.parameters["ENCUT"] == "???"
    assert any("unknown XC tag 'UNKNOWN'" in issue for issue in outcar.issues)
    assert any("ENCUT parameter value" in issue for issue in outcar.issues)


def test_outcar_compressed_and_registered(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR.bz2"
    path.write_bytes(bz2.compress(OUTCAR.encode("utf-8")))
    payload = httk.core.load(str(path), raw=True)
    assert payload["format"] == "vasp-outcar"
    assert isinstance(payload["outcar"], OutcarFile)
    assert payload["outcar"].final_energies.final is True
    assert httk.core.has_reader_for("OUTCAR")
    assert httk.core.has_reader_for("example.outcar")


def test_outcar_reader_rejects_streams() -> None:
    with pytest.raises(TypeError, match="filesystem filename"):
        read_outcar(io.StringIO(OUTCAR))


def test_outcar_prologue_is_bounded_and_construction_is_lazy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(OUTCAR + "tail\n" * 100, encoding="utf-8")
    import httk.io.vasp._text as text_module

    real_view = text_module.TextstreamFileView
    views: list[object] = []

    class CountingView:
        def __init__(self, source: object) -> None:
            self.inner = real_view(source)
            self.count = 0
            views.append(self)

        def __iter__(self) -> "CountingView":
            return self

        def __next__(self) -> str:
            self.count += 1
            return next(self.inner)

        def close(self) -> None:
            self.inner.close()

    monkeypatch.setattr(text_module, "TextstreamFileView", CountingView)
    outcar = OutcarFile(path)
    assert views == []
    assert outcar.version_string == "vasp.5.2.12"
    assert len(views) == 1
    assert views[0].count < len((OUTCAR + "tail\n" * 100).splitlines())  # type: ignore[attr-defined]


def test_outcar_closed_and_context_manager(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(OUTCAR, encoding="utf-8")
    with OutcarFile(path) as outcar:
        assert outcar.version_string == "vasp.5.2.12"
    assert outcar.closed
    with pytest.raises(ValueError, match="closed OUTCAR"):
        _ = outcar.version_string

    outcar = OutcarFile(path)
    outcar.close()
    assert outcar.closed
    with pytest.raises(ValueError, match="closed OUTCAR"):
        _ = outcar.final_energies


def test_outcar_frames_stresses_and_conversions(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(FRAME_OUTCAR, encoding="utf-8")
    outcar = OutcarFile(path)
    frames = tuple(outcar.frames())
    assert len(frames) == 3
    assert frames[0].index == 0
    assert frames[0].cell == (
        ("1.000000000", "0.000000000", "0.000000000"),
        ("0.000000000", "1.000000000", "0.000000000"),
        ("0.000000000", "0.000000000", "1.000000000"),
    )
    assert frames[1].positions == (("0.20000", "0.000000", "0.000000"), ("0.000000", "0.20000", "0.000000"))
    assert frames[2].forces == (("-0.300000", "0.000000", "0.000000"), ("0.000000", "-0.300000", "0.000000"))
    assert frames[0].free_energy == "-1.100"
    assert frames[0].stress_kbar == ("1.0", "2.0", "3.0", "4.0", "5.0", "6.0")
    assert frames[0].stress_gpa_voigt() == (-0.1, -0.2, -0.30000000000000004, -0.5, -0.6000000000000001, -0.4)
    assert outcar.nframes == 3
    assert outcar.last_frame == frames[-1]
    assert outcar.frame(1) == frames[1]
    assert outcar.stresses() == tuple(frame.stress_kbar for frame in frames)


def test_outcar_markerless_frame_stress_boundary(tmp_path: Path) -> None:
    first = frame_text(1, "-1.100", ("1.0", "2.0", "3.0", "4.0", "5.0", "6.0"))
    second = frame_text(2, "-2.200", ("2.0", "3.0", "4.0", "5.0", "6.0", "7.0")).split("\n", 1)[1]
    path = tmp_path / "markerless.outcar"
    path.write_text(FRAME_HEADER + first + second + " General timing and accounting informations\n", encoding="utf-8")
    frames = tuple(OutcarFile(path).frames())
    assert [(frame.free_energy, frame.stress_kbar) for frame in frames] == [
        ("-1.100", ("1.0", "2.0", "3.0", "4.0", "5.0", "6.0")),
        ("-2.200", ("2.0", "3.0", "4.0", "5.0", "6.0", "7.0")),
    ]


def test_outcar_replaces_rejected_builder_at_section_boundary(tmp_path: Path) -> None:
    malformed = frame_text(1, "-1.100", ("1", "2", "3", "4", "5", "6")).split("  FREE ENERGIE", 1)[0]
    second = frame_text(2, "-2.200", ("2", "3", "4", "5", "6", "7"))
    position_start = second.index("  POSITION")
    stress_start = second.index("  in kB")
    cell_start = second.index("       direct lattice vectors")
    drift_start = second.index("    total drift")
    energy_start = second.index("  FREE ENERGIE")
    markerless_second = (
        second[position_start:drift_start]
        + second[stress_start:cell_start]
        + second[drift_start:energy_start]
        + second[energy_start:]
    )
    path = tmp_path / "rejected-section.outcar"
    path.write_text(
        FRAME_HEADER + malformed + markerless_second + " General timing and accounting informations\n",
        encoding="utf-8",
    )
    frames = tuple(OutcarFile(path).frames())
    assert len(frames) == 1
    assert frames[0].free_energy == "-2.200"
    assert frames[0].stress_kbar == ("2", "3", "4", "5", "6", "7")
    assert frames[0].cell is None


def test_outcar_ions_per_type_and_nions_row_validation(tmp_path: Path) -> None:
    path = tmp_path / "ions.outcar"
    path.write_text(
        NIONS_FRAME_HEADER + frame_text(1, "-1.100", ("1", "2", "3", "4", "5", "6")) + " General timing\n",
        encoding="utf-8",
    )
    outcar = OutcarFile(path)
    assert outcar.ions_per_type == (1, 1)
    assert outcar.nframes == 1

    row = "     0.000000  0.20000  0.000000        0.000000 -0.200000  0.000000\n"
    missing_middle = frame_text(2, "-2.200", ("2", "3", "4", "5", "6", "7")).replace(row, "")
    middle_path = tmp_path / "missing-middle.outcar"
    middle_path.write_text(
        NIONS_FRAME_HEADER
        + frame_text(1, "-1.100", ("1", "2", "3", "4", "5", "6"))
        + missing_middle
        + frame_text(3, "-3.300", ("3", "4", "5", "6", "7", "8"))
        + " General timing\n",
        encoding="utf-8",
    )
    middle = OutcarFile(middle_path)
    assert middle.nframes == 2
    assert any("frame 1" in issue and "1, expected 2" in issue for issue in middle.issues)

    missing_forces = frame_text(1, "-1.100", ("1", "2", "3", "4", "5", "6")).replace(
        "     0.000000  0.10000  0.000000        0.000000 -0.100000  0.000000\n",
        "     0.000000  0.10000  0.000000\n",
    )
    final_path = tmp_path / "missing-final.outcar"
    final_path.write_text(NIONS_FRAME_HEADER + missing_forces + " General timing\n", encoding="utf-8")
    final = OutcarFile(final_path)
    assert final.nframes == 0
    assert any("frame 0" in issue and "1, expected 2" in issue for issue in final.issues)


def test_outcar_md_temperature_and_kpoint_elision(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(
        FRAME_HEADER
        + frame_text(1, "-4.400", ("1", "2", "3", "4", "5", "6"), temperature="299.68")
        + " General timing and accounting informations\n",
        encoding="utf-8",
    )
    frame = next(OutcarFile(path).frames())
    assert frame.temperature == "299.68"


def test_outcar_truncated_frame_is_dropped_and_reported(tmp_path: Path) -> None:
    complete = frame_text(1, "-1.100", ("1", "2", "3", "4", "5", "6"))
    truncated = frame_text(2, "-2.200", ("2", "3", "4", "5", "6", "7")).split("  FREE ENERGIE")[0]
    path = tmp_path / "OUTCAR"
    path.write_text(FRAME_HEADER + complete + truncated, encoding="utf-8")
    outcar = OutcarFile(path)
    assert len(tuple(outcar.frames())) == 1
    assert outcar.nframes == 1
    assert sum("truncated ionic frame" in issue for issue in outcar.issues) == 1


ELASTIC = """ TOTAL ELASTIC MODULI (kBar)
             xx          yy          zz          xy          yz          zx
 xx       1.0 2.0 3.0 4.0 5.0 6.0
 yy       7.0 8.0 9.0 10.0 11.0 12.0
 zz       13.0 14.0 15.0 16.0 17.0 18.0
 xy       19.0 20.0 21.0 22.0 23.0 24.0
 yz       25.0 26.0 27.0 28.0 29.0 30.0
 zx       31.0 32.0 33.0 34.0 35.0 36.0
 SYMMETRIZED ELASTIC MODULI (kBar)
 xx       101 102 103 104 105 106
 yy       107 108 109 110 111 112
 zz       113 114 115 116 117 118
 xy       119 120 121 122 123 124
 yz       125 126 127 128 129 130
 zx       131 132 133 134 135 136
 General timing and accounting informations
"""


def test_outcar_elastic_moduli_blocks(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(FRAME_HEADER + ELASTIC, encoding="utf-8")
    blocks = OutcarFile(path).elastic_moduli
    assert blocks == (
        ElasticModuliBlock(
            " TOTAL ELASTIC MODULI (kBar)",
            (
                ("1.0", "2.0", "3.0", "4.0", "5.0", "6.0"),
                ("7.0", "8.0", "9.0", "10.0", "11.0", "12.0"),
                ("13.0", "14.0", "15.0", "16.0", "17.0", "18.0"),
                ("19.0", "20.0", "21.0", "22.0", "23.0", "24.0"),
                ("25.0", "26.0", "27.0", "28.0", "29.0", "30.0"),
                ("31.0", "32.0", "33.0", "34.0", "35.0", "36.0"),
            ),
        ),
        ElasticModuliBlock(
            " SYMMETRIZED ELASTIC MODULI (kBar)",
            (
                ("101", "102", "103", "104", "105", "106"),
                ("107", "108", "109", "110", "111", "112"),
                ("113", "114", "115", "116", "117", "118"),
                ("119", "120", "121", "122", "123", "124"),
                ("125", "126", "127", "128", "129", "130"),
                ("131", "132", "133", "134", "135", "136"),
            ),
        ),
    )


def magnetization_block(totals: tuple[str, ...], *, closed: bool = True, tot: bool = True) -> str:
    header = (
        " magnetization (x)\n\n# of ion       s       p       d       tot\n------------------------------------------\n"
    )
    rows = "".join(f"   {i:2d}        0.000   0.000   0.000   {value}\n" for i, value in enumerate(totals, 1))
    closing = "--------------------------------------------------\n" if closed else ""
    tot_line = f"tot          0.000   0.000   0.000   {totals[-1]}\n" if tot else ""
    return header + rows + closing + tot_line


def test_outcar_magnetization_last_block_totals(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(
        FRAME_HEADER
        + magnetization_block(("1.000", "-1.000"))
        + magnetization_block(("2.734", "-2.734"))
        + " General timing and accounting informations\n",
        encoding="utf-8",
    )
    outcar = OutcarFile(path)
    assert outcar.magnetization == (2.734, -2.734)
    assert outcar.issues == ()


def test_outcar_magnetization_absent_is_none(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(OUTCAR, encoding="utf-8")
    assert OutcarFile(path).magnetization is None


def test_outcar_magnetization_missing_closing_separator_is_none(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(FRAME_HEADER + magnetization_block(("2.734", "-2.734"), closed=False), encoding="utf-8")
    outcar = OutcarFile(path)
    assert outcar.magnetization is None
    assert any("malformed magnetization (x) block" in issue for issue in outcar.issues)


def test_outcar_magnetization_missing_tot_line_is_none(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(
        FRAME_HEADER
        + magnetization_block(("2.734", "-2.734"), tot=False)
        + " General timing and accounting informations\n",
        encoding="utf-8",
    )
    outcar = OutcarFile(path)
    assert outcar.magnetization is None
    assert any("malformed magnetization (x) block" in issue for issue in outcar.issues)


def test_outcar_magnetization_noncollinear_keeps_x_totals(tmp_path: Path) -> None:
    x_block = magnetization_block(("2.734", "-2.734"))
    y_block = x_block.replace("magnetization (x)", "magnetization (y)")
    z_block = x_block.replace("magnetization (x)", "magnetization (z)")
    path = tmp_path / "OUTCAR"
    path.write_text(
        FRAME_HEADER + x_block + y_block + z_block + " General timing and accounting informations\n",
        encoding="utf-8",
    )
    outcar = OutcarFile(path)
    assert outcar.magnetization == (2.734, -2.734)
    assert outcar.noncollinear_magnetization is True


def test_outcar_magnetization_collinear_is_not_flagged(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(
        FRAME_HEADER + magnetization_block(("2.734", "-2.734")) + " General timing\n",
        encoding="utf-8",
    )
    assert OutcarFile(path).noncollinear_magnetization is False


def test_outcar_noncollinear_false_without_any_block(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    path.write_text(OUTCAR, encoding="utf-8")
    outcar = OutcarFile(path)
    assert outcar.noncollinear_magnetization is False
    assert outcar.magnetization is None


def test_outcar_frames_stream_without_full_cache(tmp_path: Path) -> None:
    path = tmp_path / "OUTCAR"
    text = FRAME_HEADER + "".join(frame_text(i, f"-{i}.000", ("1", "2", "3", "4", "5", "6")) for i in range(1, 501))
    path.write_text(text, encoding="utf-8")
    outcar = OutcarFile(path)
    assert sum(1 for _ in outcar.frames()) == 500
    assert outcar._full is None
