from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

from sprite_mask import __version__
from sprite_mask.config import AlignmentRunConfig, VcfRunConfig
from sprite_mask.workflow import run_workflow


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
        threshold=args.threshold,
        out_dir=Path(args.out),
        work_dir=Path(args.work) if args.work else None,
        threads=args.threads,
        jobs=args.jobs,
        targets_bed=Path(args.targets) if args.targets else None,
        mapq=args.mapq,
        exclude_flag=args.exclude_flag,
        reference=Path(args.reference) if args.reference else None,
        strict_depth=args.strict_depth,
        keep_work=args.keep_work,
        force=args.force,
        dry_run=args.dry_run,
    )
    outputs = run_workflow(config)
    if not args.dry_run:
        print(f"Wrote {outputs.population_count_bed_gz}")
        print(f"Wrote {outputs.population_count_bed_index}")
    return 0


def _cmd_from_vcf(args: argparse.Namespace) -> int:
    config = VcfRunConfig(
        all_sites_vcf=Path(args.all_sites_vcf),
        popfile_path=Path(args.popfile),
        threshold=args.threshold,
        out_dir=Path(args.out),
        work_dir=Path(args.work) if args.work else None,
        targets_bed=Path(args.targets) if args.targets else None,
        keep_work=args.keep_work,
        force=args.force,
        dry_run=args.dry_run,
    )
    outputs = run_workflow(config)
    if not args.dry_run:
        print(f"Wrote {outputs.population_count_bed_gz}")
        print(f"Wrote {outputs.population_count_bed_index}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sprite", description="build cohort depth mask BEDs")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(metavar="COMMAND")

    _build_from_alignments_parser(subparsers)
    _build_from_vcf_parser(subparsers)

    return parser


def _add_common_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--threshold", required=True, type=int, help="minimum passing depth")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--work", help="working directory; defaults to <out>/work")
    p.add_argument("--targets", help="optional BED regions to include")
    p.add_argument("--keep-work", action="store_true", help="keep intermediate files")
    p.add_argument("--force", action="store_true", help="overwrite existing final outputs")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="validate inputs and print plan without running",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="enable debug logging")
    p.add_argument("--quiet", "-q", action="store_true", help="suppress informational logging")
    p.add_argument(
        "--debug", action="store_true", help="re-raise exceptions instead of printing them"
    )


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
    _add_common_run_args(p)
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "mosdepth threads per sample "
            "(effective parallelism = --jobs × --threads)"
        ),
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="samples to process concurrently (effective parallelism = --jobs × --threads)",
    )
    p.add_argument("--mapq", type=int, help="mosdepth mapping-quality threshold")
    p.add_argument("--exclude-flag", type=int, help="mosdepth read flag exclusion value")
    p.add_argument("--reference", help="FASTA reference for CRAM inputs")
    p.add_argument("--strict-depth", action="store_true", help="omit mosdepth --fast-mode")
    p.set_defaults(subcommand=_cmd_from_alignments)


def _build_from_vcf_parser(subparsers: argparse._SubParsersAction) -> None:  # type: ignore[type-arg]
    p = subparsers.add_parser(
        "from-vcf",
        help="build mask from a prefiltered all-sites VCF with per-sample FORMAT/DP values",
    )
    p.add_argument(
        "--all-sites-vcf",
        required=True,
        help="prefiltered all-sites VCF with per-sample FORMAT/DP values",
    )
    p.add_argument("--popfile", required=True, help="sample/population TSV (sample_id, population)")
    _add_common_run_args(p)
    p.set_defaults(subcommand=_cmd_from_vcf)


def _setup_logging(*, verbose: bool, quiet: bool) -> None:
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s", stream=sys.stderr)


def _print_subprocess_error(error: subprocess.CalledProcessError) -> None:
    print(f"sprite: command failed with exit code {error.returncode}", file=sys.stderr)
    print(" ".join(str(part) for part in error.cmd), file=sys.stderr)
    if error.stderr:
        print(str(error.stderr).rstrip(), file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
