from __future__ import annotations

import argparse
from pathlib import Path

from sprite_mask.bedtools import build_multiinter_command
from sprite_mask.cli import build_parser
from sprite_mask.config import AlignmentRunConfig, VcfRunConfig
from sprite_mask.models import Sample
from sprite_mask.mosdepth import build_mosdepth_command
from sprite_mask.workflow import _required_tools, workflow_output_paths


def test_build_mosdepth_command_default_omits_fast_mode(tmp_path: Path) -> None:
    sample = Sample("s1", "popA", tmp_path / "s1.bam")
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=30,
        out_dir=tmp_path / "out",
        threads=4,
        min_mapq=20,
        exclude_flag=1796,
        reference=tmp_path / "ref.fa",
    )

    command = build_mosdepth_command(sample, config, tmp_path / "work" / "s1.d30")

    assert command == [
        "mosdepth",
        "--threads",
        "4",
        "--no-per-base",
        "--quantize",
        "0:30:",
        "--mapq",
        "20",
        "--flag",
        "1796",
        "--fasta",
        str(tmp_path / "ref.fa"),
        str(tmp_path / "work" / "s1.d30"),
        str(tmp_path / "s1.bam"),
    ]


def test_build_mosdepth_command_fast_mode_adds_fast_mode(tmp_path: Path) -> None:
    sample = Sample("s1", "popA", tmp_path / "s1.bam")
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=30,
        out_dir=tmp_path / "out",
        fast_mode=True,
    )

    command = build_mosdepth_command(sample, config, tmp_path / "work" / "s1.d30")

    assert "--fast-mode" in command


def test_build_mosdepth_command_max_dp_uses_three_bin_quantize(tmp_path: Path) -> None:
    sample = Sample("s1", "popA", tmp_path / "s1.bam")
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=30,
        max_dp=500,
        out_dir=tmp_path / "out",
    )

    command = build_mosdepth_command(sample, config, tmp_path / "work" / "s1.d30")

    quantize_idx = command.index("--quantize")
    assert command[quantize_idx + 1] == "0:30:500:"


def test_build_multiinter_command_uses_names_before_inputs(tmp_path: Path) -> None:
    command = build_multiinter_command(
        [tmp_path / "s1.bed", tmp_path / "s2.bed"],
        ["s1", "s2"],
    )

    assert command == [
        "bedtools",
        "multiinter",
        "-header",
        "-names",
        "s1",
        "s2",
        "-i",
        str(tmp_path / "s1.bed"),
        str(tmp_path / "s2.bed"),
    ]


def test_workflow_output_paths(tmp_path: Path) -> None:
    outputs = workflow_output_paths(tmp_path / "results", 30)

    assert outputs.population_count_bed_gz == (
        tmp_path / "results" / "sprite.bed.gz"
    )
    assert outputs.population_count_bed_index == (
        tmp_path / "results" / "sprite.bed.gz.tbi"
    )


def test_workflow_output_paths_accepts_prefix(tmp_path: Path) -> None:
    outputs = workflow_output_paths(tmp_path / "results", 30, "custom")

    assert outputs.population_count_bed_gz == tmp_path / "results" / "custom.bed.gz"
    assert outputs.population_count_bed_index == tmp_path / "results" / "custom.bed.gz.tbi"


def test_required_tools_for_vcf_mode_omits_mosdepth(tmp_path: Path) -> None:
    config = VcfRunConfig(
        all_sites_vcf=tmp_path / "all_sites.vcf.gz",
        popfile_path=tmp_path / "popfile.tsv",
        min_dp=30,
        out_dir=tmp_path / "out",
    )

    assert _required_tools(config) == ("bgzip", "tabix")


def test_cli_parser_exposes_only_run_subcommands() -> None:
    parser = build_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, argparse._SubParsersAction)
    )

    assert set(subparsers.choices) == {"from-alignments", "from-vcf"}


def test_cli_parser_from_alignments_subcommand(tmp_path: Path) -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--threads",
            "4",
            "--jobs",
            "2",
        ]
    )

    assert args.samples == str(tmp_path / "samples.tsv")
    assert args.min_dp == 30
    assert args.out == str(tmp_path / "out")
    assert args.threads == 4
    assert args.jobs == 2


def test_cli_parser_from_vcf_subcommand(tmp_path: Path) -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "from-vcf",
            "--all-sites-vcf",
            str(tmp_path / "all_sites.vcf.gz"),
            "--popfile",
            str(tmp_path / "popfile.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--snps-only",
        ]
    )

    assert args.all_sites_vcf == str(tmp_path / "all_sites.vcf.gz")
    assert args.popfile == str(tmp_path / "popfile.tsv")
    assert args.min_dp == 30
    assert args.snps_only is True
