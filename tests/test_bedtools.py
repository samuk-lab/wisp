from __future__ import annotations

import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

import pytest

from sprite_mask.bedtools import (
    build_multiinter_command,
    intersect_sort_merge,
    run_multiinter,
    sort_and_merge_bed,
    subtract_sort_merge,
    write_single_input_multiinter,
)


def test_sort_and_merge_bed_builds_bedtools_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[list[list[str]], Path]] = []

    def fake_run_pipeline(commands: Sequence[Sequence[str]], out_path: Path) -> None:
        calls.append(([list(command) for command in commands], out_path))
        out_path.write_text("merged\n")

    monkeypatch.setattr("sprite_mask.bedtools.run_pipeline", fake_run_pipeline)

    in_bed = tmp_path / "in.bed"
    out_bed = tmp_path / "out.bed"

    assert sort_and_merge_bed(in_bed, out_bed) == out_bed
    assert out_bed.read_text() == "merged\n"
    assert calls == [
        (
            [
                ["bedtools", "sort", "-i", str(in_bed)],
                ["bedtools", "merge", "-i", "-"],
            ],
            out_bed,
        )
    ]


def test_intersect_sort_merge_builds_bedtools_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[list[list[str]], Path]] = []

    def fake_run_pipeline(commands: Sequence[Sequence[str]], out_path: Path) -> None:
        calls.append(([list(command) for command in commands], out_path))

    monkeypatch.setattr("sprite_mask.bedtools.run_pipeline", fake_run_pipeline)

    a_bed = tmp_path / "a.bed"
    b_bed = tmp_path / "b.bed"
    out_bed = tmp_path / "out.bed"

    assert intersect_sort_merge(a_bed, b_bed, out_bed) == out_bed
    assert calls == [
        (
            [
                ["bedtools", "intersect", "-a", str(a_bed), "-b", str(b_bed)],
                ["bedtools", "sort", "-i", "-"],
                ["bedtools", "merge", "-i", "-"],
            ],
            out_bed,
        )
    ]


def test_subtract_sort_merge_builds_bedtools_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[list[list[str]], Path]] = []

    def fake_run_pipeline(commands: Sequence[Sequence[str]], out_path: Path) -> None:
        calls.append(([list(command) for command in commands], out_path))

    monkeypatch.setattr("sprite_mask.bedtools.run_pipeline", fake_run_pipeline)

    a_bed = tmp_path / "a.bed"
    b_bed = tmp_path / "b.bed"
    out_bed = tmp_path / "out.bed"

    assert subtract_sort_merge(a_bed, b_bed, out_bed) == out_bed
    assert calls == [
        (
            [
                ["bedtools", "subtract", "-a", str(a_bed), "-b", str(b_bed)],
                ["bedtools", "sort", "-i", "-"],
                ["bedtools", "merge", "-i", "-"],
            ],
            out_bed,
        )
    ]


def test_run_multiinter_invokes_bedtools_and_writes_stdout(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured_command: list[str] | None = None

    def fake_run(
        command: list[str],
        *,
        check: bool,
        stdout: TextIO,
        text: bool,
    ) -> subprocess.CompletedProcess[str]:
        nonlocal captured_command
        captured_command = command
        assert check is True
        assert text is True
        stdout.write("chrom\tstart\tend\tnum\tlist\ts1\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("sprite_mask.bedtools.subprocess.run", fake_run)

    out_tsv = tmp_path / "nested" / "multiinter.tsv"

    assert run_multiinter([tmp_path / "s1.bed"], ["s1"], out_tsv) == out_tsv
    assert captured_command == build_multiinter_command([tmp_path / "s1.bed"], ["s1"])
    assert out_tsv.read_text() == "chrom\tstart\tend\tnum\tlist\ts1\n"


def test_write_single_input_multiinter_skips_bed_headers_and_comments(tmp_path: Path) -> None:
    pass_bed = tmp_path / "pass.bed"
    pass_bed.write_text(
        "# comment\n"
        "chrom\tstart\tend\n"
        "chr1\t0\t10\n"
        "\n"
        "chr1\t20\t30\n"
    )
    out_tsv = tmp_path / "multiinter.tsv"

    write_single_input_multiinter(pass_bed, "s1", out_tsv)

    assert out_tsv.read_text() == (
        "chrom\tstart\tend\tnum\tlist\ts1\n"
        "chr1\t0\t10\t1\ts1\t1\n"
        "chr1\t20\t30\t1\ts1\t1\n"
    )


def test_build_multiinter_command_rejects_invalid_inputs(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="at least one"):
        build_multiinter_command([], [])

    with pytest.raises(ValueError, match="must match"):
        build_multiinter_command([tmp_path / "s1.bed"], ["s1", "s2"])
