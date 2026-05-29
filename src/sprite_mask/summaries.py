from __future__ import annotations

import gzip
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def summarize_population_count_bed(population_count_bed: Path, output_tsv: Path) -> Path:
    with _open_text(population_count_bed) as source:
        header = _read_header(source, population_count_bed)
        populations = header[3:]
        counts: dict[str, Counter[int]] = {population: Counter() for population in populations}
        for fields in _iter_interval_fields(source):
            start = int(fields[1])
            end = int(fields[2])
            length = end - start
            for population, value in zip(populations, fields[3:], strict=True):
                counts[population][int(value)] += length

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w") as out:
        out.write("population\tpassing_samples\tsites\n")
        for population in populations:
            for passing_samples in sorted(counts[population]):
                out.write(f"{population}\t{passing_samples}\t{counts[population][passing_samples]}\n")
    return output_tsv


def _open_text(path: Path) -> TextIO:
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(path, "rt")
    return path.open()


def _iter_interval_fields(source: TextIO) -> Iterator[list[str]]:
    for line in source:
        fields = line.rstrip("\n").split("\t")
        if not fields or fields == [""]:
            continue
        if fields[0].startswith("#"):
            continue
        if _looks_like_header(fields):
            continue
        yield fields


def _read_header(source: TextIO, path: Path) -> list[str]:
    for line in source:
        fields = line.rstrip("\n").split("\t")
        if not fields or fields == [""]:
            continue
        if fields[0].startswith("#"):
            header_fields = [fields[0].lstrip("#"), *fields[1:]]
            if _looks_like_header(header_fields):
                return header_fields
            continue
        if not _looks_like_header(fields):
            raise ValueError(f"{path} must have a header")
        return fields
    raise ValueError(f"{path} is empty")


def _looks_like_header(fields: list[str]) -> bool:
    return len(fields) >= 3 and fields[0] == "chrom" and fields[1] == "start" and fields[2] == "end"
