from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from sprite_mask import __version__
from sprite_mask.config import AlignmentRunConfig, VcfRunConfig
from sprite_mask.workflow import run_workflow

HELP_BANNER = """\
                     █▄
             ▄    ▀▀▄██▄
 ▄██▀█ ████▄ ████▄██ ██ ▄█▀█▄
 ▀███▄ ██ ██ ██   ██ ██ ██▄█▀
█▄▄██▀▄████▀▄█▀  ▄██▄██▄▀█▄▄▄
       ██
       ▀
"""


class SpriteArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args: Any, show_banner: bool = False, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._show_banner = show_banner

    def format_help(self) -> str:
        help_text = super().format_help()
        if not self._show_banner:
            return help_text
        return f"{HELP_BANNER}\nsprite {__version__}\n\n{help_text}"


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "subcommand"):
        parser.print_help(sys.stderr)
        return 1

    _setup_logging(verbose=getattr(args, "verbose", False), quiet=getattr(args, "quiet", False))

    try:
        return args.subcommand(args)
    except subprocess.CalledProcessError as error:
        _print_subprocess_error(error)
        return 1
    except Exception as error:
        if getattr(args, "debug", False):
            raise
        print(f"sprite: error: {error}", file=sys.stderr)
        return 1


def _cmd_from_alignments(args: argparse.Namespace) -> int:
    config = AlignmentRunConfig(
        samples_path=Path(args.samples),
        min_dp=args.min_dp,
        out_dir=Path(args.out),
        output_prefix=args.output_prefix,
        work_dir=Path(args.work) if args.work else None,
        threads=args.threads,
        jobs=args.jobs,
        mask_bed=Path(args.mask) if args.mask else None,
        variants_vcf=Path(args.variants_vcf) if args.variants_vcf else None,
        min_mapq=args.min_mapq,
        max_dp=args.max_dp,
        exclude_flag=args.exclude_flag,
        reference=Path(args.reference) if args.reference else None,
        fast_mode=args.fast_mode,
        keep_work=args.keep_work,
        force=args.force,
        dry_run=args.dry_run,
    )
    run_workflow(config)
    return 0


def _cmd_from_vcf(args: argparse.Namespace) -> int:
    config = VcfRunConfig(
        all_sites_vcf=Path(args.all_sites_vcf),
        popfile_path=Path(args.popfile),
        min_dp=args.min_dp,
        out_dir=Path(args.out),
        output_prefix=args.output_prefix,
        work_dir=Path(args.work) if args.work else None,
        mask_bed=Path(args.mask) if args.mask else None,
        keep_work=args.keep_work,
        force=args.force,
        dry_run=args.dry_run,
        snps_only=args.snps_only,
    )
    run_workflow(config)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = SpriteArgumentParser(
        prog="sprite",
        description="build population count masks",
        show_banner=True,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(
        metavar="COMMAND",
        parser_class=argparse.ArgumentParser,
    )

    _build_from_alignments_parser(subparsers)
    _build_from_vcf_parser(subparsers)

    return parser


def _add_common_run_args(p: argparse.ArgumentParser, *, min_dp_required: bool = True) -> None:
    p.add_argument(
        "--min-dp",
        required=min_dp_required,
        type=int,
        help="minimum depth to pass a site",
    )
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument(
        "--output-prefix",
        default="sprite",
        help="output file prefix within --out; .bed.gz is appended",
    )
    p.add_argument("--work", help="working directory; defaults to <out>/work")
    p.add_argument("--mask", help="restrict output to sites within this BED")
    p.add_argument("--keep-work", action="store_true", help="keep intermediate files")
    p.add_argument("--force", action="store_true", help="overwrite existing final outputs")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="validate inputs and print plan without running",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", "-q", action="store_true", help="suppress informational logging")
    p.add_argument("--debug", action="store_true", help="show full error tracebacks")


def _build_from_alignments_parser(subparsers: argparse._SubParsersAction) -> None:  # type: ignore[type-arg]
    p = subparsers.add_parser(
        "from-alignments",
        help="build mask from BAM/CRAM files via mosdepth",
    )
    p.add_argument(
        "--samples",
        required=True,
        help="sample metadata TSV (sample_id, population, alignment)",
    )
    _add_common_run_args(p, min_dp_required=False)
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="mosdepth threads per sample",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="samples to process concurrently; total parallelism = --jobs × --threads",
    )
    p.add_argument(
        "--variants-vcf",
        help=(
            "variants-only VCF used to estimate omitted depth/MAPQ thresholds "
            "and mask non-SNP variant spans"
        ),
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        help="minimum read mapping quality; defaults from --variants-vcf when available",
    )
    p.add_argument(
        "--max-dp",
        type=int,
        help="maximum depth to pass a site; defaults from --variants-vcf when available",
    )
    p.add_argument("--exclude-flag", type=int, help="SAM FLAG bits to exclude reads")
    p.add_argument("--reference", help="FASTA reference for CRAM inputs")
    p.add_argument(
        "--fast-mode",
        action="store_true",
        help=(
            "opt into mosdepth --fast-mode; strict per-base depth counting "
            "is the default"
        ),
    )
    p.set_defaults(subcommand=_cmd_from_alignments)


def _build_from_vcf_parser(subparsers: argparse._SubParsersAction) -> None:  # type: ignore[type-arg]
    p = subparsers.add_parser(
        "from-vcf",
        help="build mask from a prefiltered all-sites VCF with per-sample FORMAT/DP values",
    )
    p.add_argument(
        "--all-sites-vcf",
        required=True,
        help="all-sites VCF with per-sample FORMAT/DP",
    )
    p.add_argument("--popfile", required=True, help="sample/population TSV (sample_id, population)")
    _add_common_run_args(p)
    p.add_argument(
        "--snps-only",
        action="store_true",
        help="exclude indel sites; retain invariant sites",
    )
    p.set_defaults(subcommand=_cmd_from_vcf)


def _setup_logging(*, verbose: bool, quiet: bool) -> None:
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [sprite] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )


def _print_subprocess_error(error: subprocess.CalledProcessError) -> None:
    print(f"sprite: command failed with exit code {error.returncode}", file=sys.stderr)
    print(" ".join(str(part) for part in error.cmd), file=sys.stderr)
    if error.stderr:
        print(str(error.stderr).rstrip(), file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
