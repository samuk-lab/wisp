from __future__ import annotations

import gzip
import logging
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass
from math import floor
from pathlib import Path
from typing import TextIO

from sprite_mask.bedio import iter_bed3
from sprite_mask.collapse import write_quantized_bed_header
from sprite_mask.models import Sample
from sprite_mask.samples import populations_in_order

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class VariantThresholdEstimates:
    min_dp: int | None
    max_dp: int | None
    min_mapq: int | None
    depth_value_count: int
    mapq_value_count: int


def build_population_counts_from_all_sites_vcf(
    samples: list[Sample],
    all_sites_vcf: Path,
    output_bed: Path,
    *,
    threshold: int,
    targets_bed: Path | None = None,
    depth_field: str = "DP",
    snps_only: bool = False,
    metadata: dict[str, object] | None = None,
) -> Path:
    """Build a population-count BED directly from a prefiltered all-sites VCF."""
    target_index = _TargetIndex.from_bed(targets_bed) if targets_bed is not None else None
    populations = populations_in_order(samples)
    population_index = {population: index for index, population in enumerate(populations)}
    population_sample_counts = Counter(sample.population for sample in samples)
    sample_population_indexes = [population_index[sample.population] for sample in samples]

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(all_sites_vcf) as source, output_bed.open("w") as out:
        vcf_samples, next_line_number = _read_vcf_samples(source, all_sites_vcf)
        selected_sample_columns = _selected_sample_columns(samples, vcf_samples, all_sites_vcf)

        write_quantized_bed_header(
            out,
            columns=["chrom", "start", "end", *populations],
            metadata={
                "format": "population_count_quantized_bed",
                "source_mode": "all_sites_vcf",
                "populations": populations,
                "population_columns": [
                    {"column_number": index + 4, "name": population}
                    for index, population in enumerate(populations)
                ],
                "population_sample_counts": {
                    population: population_sample_counts[population]
                    for population in populations
                },
                "snps_only": snps_only,
                **(metadata or {}),
            },
        )

        current_coord: tuple[str, int] | None = None
        current_passes: list[bool] | None = None
        current_interval: tuple[str, int, int, tuple[int, ...]] | None = None
        last_chrom: str | None = None
        last_pos = 0
        closed_chroms: set[str] = set()

        for line_number, line in enumerate(source, start=next_line_number):
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            chrom, pos = _parse_vcf_coordinate(fields, all_sites_vcf, line_number)
            coord = (chrom, pos)
            last_chrom, last_pos = _validate_record_order(
                chrom,
                pos,
                last_chrom,
                last_pos,
                closed_chroms,
                all_sites_vcf,
                line_number,
            )

            if current_coord is not None and coord != current_coord:
                current_interval = _append_site_counts(
                    out,
                    current_interval,
                    current_coord,
                    current_passes,
                    sample_population_indexes,
                    len(populations),
                    target_index,
                )
                current_coord = None
                current_passes = None

            if snps_only and not _is_snp_or_invariant_record(fields):
                continue

            if current_passes is None:
                current_coord = coord
                current_passes = [False] * len(samples)

            depth_index = _depth_format_index(fields, depth_field, all_sites_vcf, line_number)
            _update_sample_passes(
                current_passes,
                fields,
                selected_sample_columns,
                depth_index,
                threshold,
                all_sites_vcf,
                line_number,
            )

        current_interval = _append_site_counts(
            out,
            current_interval,
            current_coord,
            current_passes,
            sample_population_indexes,
            len(populations),
            target_index,
        )
        _flush_population_count(out, current_interval)

    return output_bed


def validate_vcf_sample_names(samples: list[Sample], all_sites_vcf: Path) -> None:
    with _open_text(all_sites_vcf) as source:
        vcf_samples, _next_line_number = _read_vcf_samples(source, all_sites_vcf)
    _selected_sample_columns(samples, vcf_samples, all_sites_vcf, warn_extra=True)


def estimate_alignment_thresholds_from_variants_vcf(
    variants_vcf: Path,
    samples: list[Sample] | None = None,
) -> VariantThresholdEstimates:
    min_dp: int | None = None
    max_dp: int | None = None
    min_mq: float | None = None
    depth_value_count = 0
    mapq_value_count = 0

    with _open_text(variants_vcf) as source:
        vcf_samples, next_line_number, _header_line_number = _read_vcf_samples_optional(
            source,
            variants_vcf,
        )
        selected_sample_columns: list[int] | None = None
        if samples is not None and vcf_samples:
            selected_sample_columns = _selected_sample_columns(samples, vcf_samples, variants_vcf)

        for line_number, line in enumerate(source, start=next_line_number):
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            _parse_vcf_fixed_coordinate(fields, variants_vcf, line_number)
            if not _has_variant_alt(fields):
                continue

            info_values = _parse_info_values(fields[7])
            for mapq in _numeric_info_values(info_values, "MQ", variants_vcf, line_number):
                min_mq = mapq if min_mq is None else min(min_mq, mapq)
                mapq_value_count += 1

            depth_index = _optional_depth_format_index(fields)
            if depth_index is None:
                continue

            sample_columns = (
                selected_sample_columns
                if selected_sample_columns is not None
                else list(range(max(0, len(fields) - 9)))
            )
            for sample_column in sample_columns:
                field_index = 9 + sample_column
                if field_index >= len(fields):
                    raise ValueError(
                        f"{variants_vcf}:{line_number} has fewer sample columns than the header"
                    )
                depth = _sample_integer_format_value(
                    fields[field_index],
                    depth_index,
                    "DP",
                    variants_vcf,
                    line_number,
                )
                if depth is None or depth <= 0:
                    continue
                min_dp = depth if min_dp is None else min(min_dp, depth)
                max_dp = depth if max_dp is None else max(max_dp, depth)
                depth_value_count += 1

    return VariantThresholdEstimates(
        min_dp=min_dp,
        max_dp=max_dp,
        min_mapq=floor(min_mq) if min_mq is not None else None,
        depth_value_count=depth_value_count,
        mapq_value_count=mapq_value_count,
    )


def write_variant_exclusion_bed(variants_vcf: Path, output_bed: Path) -> int:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    interval_count = 0

    with _open_text(variants_vcf) as source, output_bed.open("w") as out:
        _vcf_samples, next_line_number, _header_line_number = _read_vcf_samples_optional(
            source,
            variants_vcf,
        )
        for line_number, line in enumerate(source, start=next_line_number):
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            interval = _variant_exclusion_interval(fields, variants_vcf, line_number)
            if interval is None:
                continue
            chrom, start, end = interval
            out.write(f"{chrom}\t{start}\t{end}\n")
            interval_count += 1

    return interval_count


@dataclass(frozen=True)
class _TargetIndex:
    starts_by_chrom: dict[str, list[int]]
    ends_by_chrom: dict[str, list[int]]

    @classmethod
    def from_bed(cls, targets_bed: Path) -> _TargetIndex:
        intervals_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in iter_bed3(targets_bed):
            intervals_by_chrom.setdefault(chrom, []).append((start, end))

        starts_by_chrom: dict[str, list[int]] = {}
        ends_by_chrom: dict[str, list[int]] = {}
        for chrom, intervals in intervals_by_chrom.items():
            merged = _merge_intervals(sorted(intervals))
            starts_by_chrom[chrom] = [start for start, _end in merged]
            ends_by_chrom[chrom] = [end for _start, end in merged]
        return cls(starts_by_chrom=starts_by_chrom, ends_by_chrom=ends_by_chrom)

    def contains(self, chrom: str, start: int, end: int) -> bool:
        starts = self.starts_by_chrom.get(chrom)
        if starts is None:
            return False
        index = bisect_right(starts, start) - 1
        if index < 0:
            return False
        return end <= self.ends_by_chrom[chrom][index]


def _open_text(path: Path) -> TextIO:
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(path, "rt")
    return path.open()


def _read_vcf_samples(source: TextIO, path: Path) -> tuple[list[str], int]:
    sample_names, next_line_number, header_line_number = _read_vcf_samples_optional(source, path)
    if not sample_names:
        raise ValueError(f"{path}:{header_line_number} does not contain VCF sample columns")
    return sample_names, next_line_number


def _read_vcf_samples_optional(source: TextIO, path: Path) -> tuple[list[str], int, int]:
    line_number = 0
    for line_number, line in enumerate(source, start=1):
        if line.startswith("##"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields and fields[0] == "#CHROM":
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number} must have VCF fixed fields")
            sample_names = fields[9:] if len(fields) > 9 else []
            duplicates = sorted(
                sample for sample, count in Counter(sample_names).items() if count > 1
            )
            if duplicates:
                raise ValueError(
                    f"{path}:{line_number} contains duplicate VCF sample(s): "
                    + ", ".join(duplicates)
                )
            return sample_names, line_number + 1, line_number
        if line.startswith("#"):
            continue
        raise ValueError(f"{path}:{line_number} appears before the #CHROM header")
    raise ValueError(f"{path} does not contain a #CHROM header")


def _selected_sample_columns(
    samples: list[Sample],
    vcf_samples: list[str],
    path: Path,
    *,
    warn_extra: bool = False,
) -> list[int]:
    vcf_sample_columns = {sample: index for index, sample in enumerate(vcf_samples)}
    requested_sample_ids = {sample.sample_id for sample in samples}
    missing_samples = [
        sample.sample_id for sample in samples if sample.sample_id not in vcf_sample_columns
    ]
    extra_vcf_samples = [sample for sample in vcf_samples if sample not in requested_sample_ids]
    if missing_samples:
        details = "popfile sample(s) absent from VCF: " + ", ".join(missing_samples)
        raise ValueError(f"{path} sample mismatch with popfile: {details}")
    if extra_vcf_samples and warn_extra:
        logger.warning(
            "%s contains VCF sample(s) absent from popfile; these samples will be ignored: %s",
            path,
            ", ".join(extra_vcf_samples),
        )
    return [vcf_sample_columns[sample.sample_id] for sample in samples]


def _parse_vcf_coordinate(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < 10:
        raise ValueError(f"{path}:{line_number} must have VCF fixed fields and sample columns")
    return _parse_vcf_fixed_coordinate(fields, path, line_number)


def _parse_vcf_fixed_coordinate(
    fields: list[str],
    path: Path,
    line_number: int,
) -> tuple[str, int]:
    if len(fields) < 8:
        raise ValueError(f"{path}:{line_number} must have VCF fixed fields")
    chrom = fields[0]
    try:
        pos = int(fields[1])
    except ValueError as error:
        raise ValueError(f"{path}:{line_number} has a non-integer POS") from error
    if pos < 1:
        raise ValueError(f"{path}:{line_number} has POS < 1")
    return chrom, pos


def _has_variant_alt(fields: list[str]) -> bool:
    alt = fields[4]
    return alt not in {"", "."}


def _variant_exclusion_interval(
    fields: list[str],
    path: Path,
    line_number: int,
) -> tuple[str, int, int] | None:
    chrom, pos = _parse_vcf_fixed_coordinate(fields, path, line_number)
    ref = fields[3]
    alt = fields[4]
    info_values = _parse_info_values(fields[7])
    if not _requires_variant_exclusion(ref, alt, info_values):
        return None

    start = pos - 1
    end = _variant_exclusion_end(start, pos, ref, info_values, path, line_number)
    if end <= start:
        end = start + 1
    return chrom, start, end


def _requires_variant_exclusion(
    ref: str,
    alt: str,
    info_values: dict[str, list[str]],
) -> bool:
    if alt in {"", "."}:
        return False
    if "SVTYPE" in info_values:
        return True
    if len(ref) != 1 or ref in {"", ".", "*"}:
        return True

    for allele in alt.split(","):
        if allele in {"", "."}:
            continue
        if _is_symbolic_or_breakend_allele(allele):
            return True
        if len(allele) != 1:
            return True
    return False


def _variant_exclusion_end(
    start: int,
    pos: int,
    ref: str,
    info_values: dict[str, list[str]],
    path: Path,
    line_number: int,
) -> int:
    info_end = _first_integer_info_value(info_values, "END", path, line_number)
    if info_end is not None:
        if info_end < pos:
            raise ValueError(f"{path}:{line_number} has INFO/END before POS")
        return info_end

    svtype = next((value.upper() for value in info_values.get("SVTYPE", []) if value), None)
    svlens = [
        abs(value)
        for value in _integer_info_values(info_values, "SVLEN", path, line_number)
        if value != 0
    ]
    if svlens and svtype != "INS":
        return start + max(svlens)

    return start + max(1, len(ref) if ref not in {"", ".", "*"} else 1)


def _is_symbolic_or_breakend_allele(allele: str) -> bool:
    return (
        allele == "*"
        or (allele.startswith("<") and allele.endswith(">"))
        or "[" in allele
        or "]" in allele
    )


def _is_snp_or_invariant_record(fields: list[str]) -> bool:
    ref = fields[3]
    alt = fields[4]
    if len(ref) != 1 or ref in {".", "*"}:
        return False
    if alt == ".":
        return True

    alt_alleles = alt.split(",")
    return all(len(allele) == 1 and allele not in {".", "*"} for allele in alt_alleles)


def _validate_record_order(
    chrom: str,
    pos: int,
    last_chrom: str | None,
    last_pos: int,
    closed_chroms: set[str],
    path: Path,
    line_number: int,
) -> tuple[str, int]:
    if last_chrom is None:
        return chrom, pos

    if chrom != last_chrom:
        closed_chroms.add(last_chrom)
        if chrom in closed_chroms:
            raise ValueError(
                f"{path}:{line_number} has chromosome {chrom!r} in multiple blocks; "
                "VCF records must be grouped by chromosome"
            )
        return chrom, pos

    if pos < last_pos:
        raise ValueError(
            f"{path}:{line_number} has POS before an earlier {chrom!r} record; "
            "VCF records must be sorted by POS within each chromosome so duplicate "
            "CHROM:POS records are contiguous"
        )

    return chrom, pos


def _depth_format_index(
    fields: list[str],
    depth_field: str,
    path: Path,
    line_number: int,
) -> int:
    format_fields = fields[8].split(":")
    try:
        return format_fields.index(depth_field)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} FORMAT does not contain {depth_field!r}"
        ) from error


def _optional_depth_format_index(fields: list[str], depth_field: str = "DP") -> int | None:
    if len(fields) < 10:
        return None
    format_text = fields[8]
    if format_text in {"", "."}:
        return None
    format_fields = format_text.split(":")
    try:
        return format_fields.index(depth_field)
    except ValueError:
        return None


def _parse_info_values(info_text: str) -> dict[str, list[str]]:
    if info_text in {"", "."}:
        return {}

    values: dict[str, list[str]] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" not in item:
            values.setdefault(item, []).append("")
            continue
        key, value = item.split("=", maxsplit=1)
        values.setdefault(key, []).extend(value.split(","))
    return values


def _numeric_info_values(
    info_values: dict[str, list[str]],
    key: str,
    path: Path,
    line_number: int,
) -> list[float]:
    values: list[float] = []
    for value in info_values.get(key, []):
        if value in {"", "."}:
            continue
        try:
            values.append(float(value))
        except ValueError as error:
            raise ValueError(f"{path}:{line_number} has non-numeric INFO/{key}") from error
    return values


def _integer_info_values(
    info_values: dict[str, list[str]],
    key: str,
    path: Path,
    line_number: int,
) -> list[int]:
    values: list[int] = []
    for value in info_values.get(key, []):
        if value in {"", "."}:
            continue
        try:
            values.append(int(value))
        except ValueError as error:
            raise ValueError(f"{path}:{line_number} has non-integer INFO/{key}") from error
    return values


def _first_integer_info_value(
    info_values: dict[str, list[str]],
    key: str,
    path: Path,
    line_number: int,
) -> int | None:
    values = _integer_info_values(info_values, key, path, line_number)
    return values[0] if values else None


def _update_sample_passes(
    passes: list[bool],
    fields: list[str],
    selected_sample_columns: list[int],
    depth_index: int,
    threshold: int,
    path: Path,
    line_number: int,
) -> None:
    for selected_index, vcf_sample_index in enumerate(selected_sample_columns):
        if passes[selected_index]:
            continue

        field_index = 9 + vcf_sample_index
        if field_index >= len(fields):
            raise ValueError(f"{path}:{line_number} has fewer sample columns than the header")
        passes[selected_index] = _sample_depth_passes(
            fields[field_index],
            depth_index,
            threshold,
            path,
            line_number,
        )


def _sample_depth_passes(
    sample_field: str,
    depth_index: int,
    threshold: int,
    path: Path,
    line_number: int,
) -> bool:
    if sample_field in {"", "."}:
        return False
    parts = sample_field.split(":")
    if depth_index >= len(parts):
        return False

    depth_text = parts[depth_index]
    if depth_text in {"", "."}:
        return False

    try:
        return int(depth_text) >= threshold
    except ValueError as error:
        raise ValueError(f"{path}:{line_number} has non-integer sample DP") from error


def _sample_integer_format_value(
    sample_field: str,
    value_index: int,
    value_name: str,
    path: Path,
    line_number: int,
) -> int | None:
    if sample_field in {"", "."}:
        return None
    parts = sample_field.split(":")
    if value_index >= len(parts):
        return None

    value_text = parts[value_index]
    if value_text in {"", "."}:
        return None

    try:
        return int(value_text)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} has non-integer sample {value_name}"
        ) from error


def _append_site_counts(
    out: TextIO,
    current_interval: tuple[str, int, int, tuple[int, ...]] | None,
    coord: tuple[str, int] | None,
    passes: list[bool] | None,
    sample_population_indexes: list[int],
    population_count: int,
    target_index: _TargetIndex | None,
) -> tuple[str, int, int, tuple[int, ...]] | None:
    if coord is None or passes is None:
        return current_interval

    chrom, pos = coord
    start = pos - 1
    end = pos
    if target_index is not None and not target_index.contains(chrom, start, end):
        _flush_population_count(out, current_interval)
        return None

    counts = [0] * population_count
    for sample_passes, population_index in zip(
        passes,
        sample_population_indexes,
        strict=True,
    ):
        if sample_passes:
            counts[population_index] += 1
    count_tuple = tuple(counts)

    if not any(count_tuple):
        _flush_population_count(out, current_interval)
        return None

    if (
        current_interval
        and current_interval[0] == chrom
        and current_interval[2] == start
        and current_interval[3] == count_tuple
    ):
        return (current_interval[0], current_interval[1], end, current_interval[3])

    _flush_population_count(out, current_interval)
    return (chrom, start, end, count_tuple)


def _flush_population_count(
    handle: TextIO,
    current: tuple[str, int, int, tuple[int, ...]] | None,
) -> None:
    if current is None:
        return
    chrom, start, end, counts = current
    handle.write("\t".join([chrom, str(start), str(end), *[str(count) for count in counts]]) + "\n")


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
            continue
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged
