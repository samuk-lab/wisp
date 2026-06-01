from __future__ import annotations

import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from dataclasses import replace
from pathlib import Path

from sprite_mask.bedio import extract_merged_pass_intervals, normalize_targets_bed
from sprite_mask.bedtools import (
    intersect_sort_merge,
    run_multiinter,
    sort_and_merge_bed,
    subtract_sort_merge,
    write_single_input_multiinter,
)
from sprite_mask.collapse import collapse_population_counts
from sprite_mask.config import AlignmentRunConfig, RunConfig, VcfRunConfig
from sprite_mask.models import MosdepthOutputs, Sample, WorkflowOutputs
from sprite_mask.mosdepth import run_mosdepth
from sprite_mask.samples import read_popfile, read_samples
from sprite_mask.validation import (
    ensure_parent_dirs,
    refuse_existing_outputs,
    require_executables,
    validate_alignment_sample_headers,
    validate_jobs,
    validate_threads,
    validate_threshold,
    validate_variants_vcf_input,
    validate_vcf_inputs,
)
from sprite_mask.vcf import (
    build_population_counts_from_all_sites_vcf,
    estimate_alignment_thresholds_from_variants_vcf,
    validate_vcf_sample_names,
    write_variant_exclusion_bed,
)

logger = logging.getLogger(__name__)


def run_workflow(config: RunConfig) -> WorkflowOutputs:
    mode = _workflow_mode(config)
    logger.info("Analysis start: validating %s workflow inputs", mode)
    alignment_metadata: dict[str, object] = {}

    if isinstance(config, VcfRunConfig):
        logger.info("Analysis validation: checking threshold and VCF input paths")
        validate_threshold(
            config.min_dp,
            targets_bed=config.mask_bed,
            all_sites_vcf=config.all_sites_vcf,
        )
        validate_vcf_inputs(config.all_sites_vcf, config.popfile_path)
        logger.info("Analysis input: reading population file %s", config.popfile_path)
        samples = read_popfile(config.popfile_path)
    else:
        logger.info("Analysis validation: checking alignment run settings")
        validate_threads(config.threads)
        validate_jobs(config.jobs)
        logger.info("Analysis input: reading sample file %s", config.samples_path)
        samples = read_samples(config.samples_path)
        config, alignment_metadata = _resolve_alignment_thresholds(config, samples)
        logger.info("Analysis validation: checking threshold")
        assert config.min_dp is not None
        validate_threshold(config.min_dp, targets_bed=config.mask_bed)

    _log_sample_summary(samples)
    required_tools = _required_tools(config)
    logger.info("Analysis validation: checking required tools: %s", ", ".join(required_tools))
    require_executables(required_tools)
    if isinstance(config, VcfRunConfig):
        logger.info("Analysis validation: checking VCF sample columns against population file")
        validate_vcf_sample_names(samples, config.all_sites_vcf)
    else:
        logger.info("Analysis validation: checking alignment headers against sample file")
        validate_alignment_sample_headers(samples)

    assert config.min_dp is not None
    outputs = workflow_output_paths(config.out_dir, config.min_dp, config.output_prefix)
    final_paths = [
        outputs.population_count_bed_gz,
        outputs.population_count_bed_index,
    ]
    logger.info("Analysis output: checking final output paths in %s", config.out_dir)
    refuse_existing_outputs(final_paths, force=config.force)

    if config.dry_run:
        logger.info(
            "Analysis dry-run summary: would process %d sample(s) via %s mode with min-dp %d",
            len(samples),
            mode,
            config.min_dp,
        )
        logger.info(
            "Analysis dry-run summary: output would be written to %s",
            outputs.population_count_bed_gz,
        )
        logger.info("Analysis dry-run complete: no files were written")
        return outputs

    logger.info("Analysis setup: preparing output and work directories")
    ensure_parent_dirs(final_paths)
    config.resolved_work_dir.mkdir(parents=True, exist_ok=True)

    generated_work_files: list[Path] = []
    if isinstance(config, VcfRunConfig):
        _build_from_all_sites_vcf(samples, config, generated_work_files)
    else:
        _build_from_alignments(
            samples,
            config,
            generated_work_files,
            threshold_metadata=alignment_metadata,
        )

    if not config.keep_work:
        logger.info(
            "Analysis cleanup: removing temporary work files from %s",
            config.resolved_work_dir,
        )
        _cleanup_work_files(generated_work_files, config.resolved_work_dir)
    else:
        logger.info("Analysis cleanup: keeping work directory %s", config.resolved_work_dir)

    logger.info("Analysis complete: wrote %s", outputs.population_count_bed_gz)
    logger.info("Analysis complete: wrote %s", outputs.population_count_bed_index)

    return outputs


def _resolve_alignment_thresholds(
    config: AlignmentRunConfig,
    samples: list[Sample],
) -> tuple[AlignmentRunConfig, dict[str, object]]:
    if config.variants_vcf is None:
        if config.min_dp is None:
            raise ValueError("--min-dp is required unless --variants-vcf is provided")
        return config, {}

    validate_variants_vcf_input(config.variants_vcf)

    missing_thresholds = (
        config.min_dp is None or config.max_dp is None or config.min_mapq is None
    )
    estimates = None
    if missing_thresholds:
        logger.info(
            "Analysis variants VCF: estimating omitted thresholds from %s",
            config.variants_vcf,
        )
        estimates = estimate_alignment_thresholds_from_variants_vcf(
            config.variants_vcf,
            samples,
        )

    min_dp = config.min_dp
    max_dp = config.max_dp
    min_mapq = config.min_mapq
    threshold_sources = {
        "min_dp": "manual" if min_dp is not None else "variants_vcf",
        "max_dp": "manual" if max_dp is not None else "variants_vcf",
        "min_mapq": "manual" if min_mapq is not None else "variants_vcf",
    }

    if estimates is not None:
        if min_dp is None:
            min_dp = estimates.min_dp
        if max_dp is None:
            max_dp = estimates.max_dp
        if min_mapq is None:
            min_mapq = estimates.min_mapq

    if min_dp is None:
        raise ValueError(
            "--min-dp was not provided and could not be estimated from --variants-vcf; "
            "provide --min-dp or use a VCF with per-sample FORMAT/DP values"
        )

    logger.info(
        "Analysis thresholds: using min-dp=%s, max-dp=%s, min-mapq=%s",
        min_dp,
        "none" if max_dp is None else max_dp,
        "none" if min_mapq is None else min_mapq,
    )

    metadata: dict[str, object] = {
        "variants_vcf": str(config.variants_vcf),
        "threshold_sources": threshold_sources,
    }
    if estimates is not None:
        metadata["variant_vcf_threshold_estimates"] = {
            "min_dp": estimates.min_dp,
            "max_dp": estimates.max_dp,
            "min_mapq": estimates.min_mapq,
            "depth_value_count": estimates.depth_value_count,
            "mapq_value_count": estimates.mapq_value_count,
        }

    return replace(config, min_dp=min_dp, max_dp=max_dp, min_mapq=min_mapq), metadata


def _build_from_all_sites_vcf(
    samples: list[Sample],
    config: VcfRunConfig,
    generated_work_files: list[Path],
) -> None:
    metadata = {
        "min_dp": config.min_dp,
        "sample_count": len(samples),
        "popfile": str(config.popfile_path),
        "all_sites_vcf": str(config.all_sites_vcf),
        "mask_bed": str(config.mask_bed) if config.mask_bed is not None else None,
        "snps_only": config.snps_only,
    }
    population_count_bed = (
        config.resolved_work_dir / f"cohort.d{config.min_dp}.population_count_quantized.bed"
    )
    generated_work_files.append(population_count_bed)
    logger.info("Analysis VCF: building population counts from %s", config.all_sites_vcf)
    if config.mask_bed is not None:
        logger.info(
            "Analysis VCF: restricting population counts to targets in %s",
            config.mask_bed,
        )
    build_population_counts_from_all_sites_vcf(
        samples,
        config.all_sites_vcf,
        population_count_bed,
        threshold=config.min_dp,
        targets_bed=config.mask_bed,
        snps_only=config.snps_only,
        metadata=metadata,
    )

    outputs = workflow_output_paths(config.out_dir, config.min_dp, config.output_prefix)
    logger.info("Analysis VCF: preparing final indexed BED")
    _sort_bgzip_tabix_bed(population_count_bed, outputs.population_count_bed_gz)


def _build_from_alignments(
    samples: list[Sample],
    config: AlignmentRunConfig,
    generated_work_files: list[Path],
    *,
    threshold_metadata: dict[str, object] | None = None,
) -> None:
    assert config.min_dp is not None
    target_bed = _prepare_targets(config, generated_work_files)
    variant_exclusion_bed = _prepare_variant_exclusions(config, generated_work_files)
    passing_beds = _make_sample_pass_beds(
        samples,
        config,
        target_bed,
        variant_exclusion_bed,
        generated_work_files,
    )

    multiinter_tsv = config.resolved_work_dir / f"cohort.d{config.min_dp}.multiinter.tsv"
    generated_work_files.append(multiinter_tsv)
    sample_names = [sample.sample_id for sample in samples]
    logger.info("Analysis alignments: combining %d sample pass BED(s)", len(passing_beds))
    if len(passing_beds) == 1:
        write_single_input_multiinter(passing_beds[0], sample_names[0], multiinter_tsv)
    else:
        run_multiinter(passing_beds, sample_names, multiinter_tsv)

    metadata = {
        "min_dp": config.min_dp,
        "sample_count": len(samples),
        "samples_path": str(config.samples_path),
        "mask_bed": str(config.mask_bed) if config.mask_bed is not None else None,
        **(threshold_metadata or {}),
    }
    population_count_bed = (
        config.resolved_work_dir / f"cohort.d{config.min_dp}.population_count_quantized.bed"
    )
    generated_work_files.append(population_count_bed)
    logger.info("Analysis alignments: collapsing sample indicators to population counts")
    collapse_population_counts(
        samples,
        multiinter_tsv,
        population_count_bed,
        metadata=metadata,
    )
    outputs = workflow_output_paths(config.out_dir, config.min_dp, config.output_prefix)
    logger.info("Analysis alignments: preparing final indexed BED")
    _sort_bgzip_tabix_bed(population_count_bed, outputs.population_count_bed_gz)


def workflow_output_paths(
    out_dir: Path,
    threshold: int,
    output_prefix: str = "sprite",
) -> WorkflowOutputs:
    if output_prefix == "":
        raise ValueError("--output-prefix cannot be empty")
    population_count_bed_gz = out_dir / f"{output_prefix}.bed.gz"
    return WorkflowOutputs(
        population_count_bed_gz=population_count_bed_gz,
        population_count_bed_index=Path(f"{population_count_bed_gz}.tbi"),
    )


def _required_tools(config: RunConfig) -> tuple[str, ...]:
    if isinstance(config, VcfRunConfig):
        return ("bgzip", "tabix")
    if config.min_dp == 0:
        return ("samtools", "bedtools", "bgzip", "tabix")
    return ("samtools", "mosdepth", "bedtools", "bgzip", "tabix")


def _prepare_targets(config: RunConfig, generated_work_files: list[Path]) -> Path | None:
    if config.mask_bed is None:
        logger.info("Analysis targets: no target BED provided")
        return None

    normalized = config.resolved_work_dir / "targets.3col.bed"
    sorted_merged = config.resolved_work_dir / "targets.3col.sorted.merged.bed"
    generated_work_files.extend([normalized, sorted_merged])
    logger.info("Analysis targets: normalizing %s", config.mask_bed)
    normalize_targets_bed(config.mask_bed, normalized)
    logger.info("Analysis targets: sorting and merging target intervals")
    sort_and_merge_bed(normalized, sorted_merged)
    return sorted_merged


def _prepare_variant_exclusions(
    config: AlignmentRunConfig,
    generated_work_files: list[Path],
) -> Path | None:
    if config.variants_vcf is None:
        return None

    raw_exclusions = config.resolved_work_dir / "variants_vcf.excluded.raw.bed"
    merged_exclusions = config.resolved_work_dir / "variants_vcf.excluded.sorted.merged.bed"
    generated_work_files.append(raw_exclusions)

    logger.info(
        "Analysis variants VCF: writing non-SNP exclusion intervals from %s",
        config.variants_vcf,
    )
    exclusion_count = write_variant_exclusion_bed(config.variants_vcf, raw_exclusions)
    if exclusion_count == 0:
        logger.info("Analysis variants VCF: no non-SNP exclusion intervals found")
        return None

    generated_work_files.append(merged_exclusions)
    logger.info(
        "Analysis variants VCF: sorting and merging %d exclusion interval(s)",
        exclusion_count,
    )
    sort_and_merge_bed(raw_exclusions, merged_exclusions)
    return merged_exclusions


def _make_sample_pass_beds(
    samples: list[Sample],
    config: AlignmentRunConfig,
    target_bed: Path | None,
    variant_exclusion_bed: Path | None,
    generated_work_files: list[Path],
) -> list[Path]:
    logger.info(
        "Analysis alignments: processing %d sample alignment(s) "
        "with %d job(s) and %d thread(s) per sample",
        len(samples),
        min(config.jobs, len(samples)),
        config.threads,
    )
    if config.jobs == 1 or len(samples) == 1:
        results = [
            _make_sample_pass_bed(sample, config, target_bed, variant_exclusion_bed)
            for sample in samples
        ]
    else:
        max_workers = min(config.jobs, len(samples))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(
                executor.map(
                    lambda sample: _make_sample_pass_bed(
                        sample,
                        config,
                        target_bed,
                        variant_exclusion_bed,
                    ),
                    samples,
                )
            )

    passing_beds: list[Path] = []
    for pass_bed, _mosdepth_output, sample_work_files in results:
        passing_beds.append(pass_bed)
        generated_work_files.extend(sample_work_files)
    return passing_beds


def _make_sample_pass_bed(
    sample: Sample,
    config: AlignmentRunConfig,
    target_bed: Path | None,
    variant_exclusion_bed: Path | None = None,
) -> tuple[Path, MosdepthOutputs | None, list[Path]]:
    sample_prefix = config.resolved_work_dir / f"{sample.sample_id}.d{config.min_dp}"
    generated_work_files: list[Path] = []

    if config.min_dp == 0:
        if target_bed is None:
            raise ValueError("--min-dp 0 requires --mask")
        logger.info(
            "Analysis sample %s: using target intervals because threshold is 0",
            sample.sample_id,
        )
        target_copy = Path(f"{sample_prefix}.pass.targets.bed")
        shutil.copyfile(target_bed, target_copy)
        generated_work_files.append(target_copy)
        final_pass_bed = _exclude_variant_regions_from_sample_pass_bed(
            sample,
            target_copy,
            variant_exclusion_bed,
            generated_work_files,
        )
        return final_pass_bed, None, generated_work_files

    logger.info("Analysis sample %s: running mosdepth", sample.sample_id)
    mosdepth_outputs = run_mosdepth(sample, config)
    generated_work_files.extend(
        [
            mosdepth_outputs.quantized_bed_gz,
            mosdepth_outputs.quantized_bed_index,
            mosdepth_outputs.summary,
            mosdepth_outputs.global_dist,
            mosdepth_outputs.stderr_log,
        ]
    )

    merged_pass_bed = Path(f"{sample_prefix}.pass.bed")
    generated_work_files.append(merged_pass_bed)
    logger.info("Analysis sample %s: extracting pass intervals", sample.sample_id)
    extract_merged_pass_intervals(mosdepth_outputs.quantized_bed_gz, merged_pass_bed)

    if target_bed is None:
        logger.info("Analysis sample %s: finished pass interval BED", sample.sample_id)
        final_pass_bed = _exclude_variant_regions_from_sample_pass_bed(
            sample,
            merged_pass_bed,
            variant_exclusion_bed,
            generated_work_files,
        )
        return final_pass_bed, mosdepth_outputs, generated_work_files

    clipped_pass_bed = Path(f"{sample_prefix}.pass.targets.bed")
    generated_work_files.append(clipped_pass_bed)
    logger.info("Analysis sample %s: clipping pass intervals to targets", sample.sample_id)
    intersect_sort_merge(merged_pass_bed, target_bed, clipped_pass_bed)
    logger.info("Analysis sample %s: finished targeted pass interval BED", sample.sample_id)
    final_pass_bed = _exclude_variant_regions_from_sample_pass_bed(
        sample,
        clipped_pass_bed,
        variant_exclusion_bed,
        generated_work_files,
    )
    return final_pass_bed, mosdepth_outputs, generated_work_files


def _exclude_variant_regions_from_sample_pass_bed(
    sample: Sample,
    pass_bed: Path,
    variant_exclusion_bed: Path | None,
    generated_work_files: list[Path],
) -> Path:
    if variant_exclusion_bed is None:
        return pass_bed

    excluded_pass_bed = Path(f"{pass_bed.with_suffix('')}.variants.bed")
    generated_work_files.append(excluded_pass_bed)
    logger.info(
        "Analysis sample %s: removing variants-only VCF exclusion intervals",
        sample.sample_id,
    )
    subtract_sort_merge(pass_bed, variant_exclusion_bed, excluded_pass_bed)
    return excluded_pass_bed


def _sort_bgzip_tabix_bed(in_bed: Path, out_bed_gz: Path) -> None:
    sorted_bed = out_bed_gz.with_suffix("")
    body_bed = sorted_bed.with_suffix(f"{sorted_bed.suffix}.body")
    sorted_body_bed = sorted_bed.with_suffix(f"{sorted_bed.suffix}.sorted_body")
    header_lines: list[str] = []

    out_bed_gz.parent.mkdir(parents=True, exist_ok=True)
    with in_bed.open() as source, body_bed.open("w") as body:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields == [""]:
                continue
            if line.startswith("#"):
                header_lines.append(line)
                continue
            if len(fields) >= 3 and fields[0] == "chrom" and fields[1] == "start":
                header_lines.append("#" + line)
                continue
            body.write(line)

    try:
        logger.info("Analysis output: sorting BED intervals for %s", out_bed_gz)
        with sorted_body_bed.open("w") as sorted_body_out:
            subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(body_bed)],
                check=True,
                stdout=sorted_body_out,
                text=True,
                env={**__import__("os").environ, "LC_ALL": "C"},
            )

        with sorted_bed.open("w") as sorted_out, sorted_body_bed.open() as sorted_body:
            sorted_out.writelines(header_lines)
            for line in sorted_body:
                sorted_out.write(line)

        logger.info("Analysis output: compressing %s with bgzip", sorted_bed)
        subprocess.run(["bgzip", "-f", str(sorted_bed)], check=True)
        logger.info("Analysis output: indexing %s with tabix", out_bed_gz)
        subprocess.run(["tabix", "-f", "-p", "bed", str(out_bed_gz)], check=True)
    finally:
        with suppress(FileNotFoundError):
            body_bed.unlink()
        with suppress(FileNotFoundError):
            sorted_body_bed.unlink()
        with suppress(FileNotFoundError):
            sorted_bed.unlink()


def _cleanup_work_files(paths: list[Path], work_dir: Path) -> None:
    for path in reversed(paths):
        try:
            path.unlink()
        except FileNotFoundError:
            continue

    with suppress(OSError):
        work_dir.rmdir()


def _workflow_mode(config: RunConfig) -> str:
    if isinstance(config, VcfRunConfig):
        return "VCF"
    return "alignment"


def _log_sample_summary(samples: list[Sample]) -> None:
    populations: list[str] = []
    seen_populations: set[str] = set()
    for sample in samples:
        if sample.population in seen_populations:
            continue
        populations.append(sample.population)
        seen_populations.add(sample.population)

    logger.info(
        "Analysis input summary: loaded %d sample(s) across %d population(s): %s",
        len(samples),
        len(populations),
        ", ".join(populations),
    )
