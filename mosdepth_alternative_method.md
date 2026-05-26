# Mosdepth alternative method

## Goal

For a fixed depth threshold, calculate how many sites in each BAM clear that threshold, then intersect the passing intervals across BAMs to count how many samples meet the threshold at each genomic site. Finally, collapse adjacent windows with the same sample count into ranges, similar to a quantized BED file.

The same workflow can also produce population-aware count beds, where each output interval has one column per population giving the number of samples from that population that pass the depth threshold.

Core tools: [`mosdepth`](https://github.com/brentp/mosdepth) for per-sample threshold pass intervals and [`bedtools multiinter`](https://bedtools.readthedocs.io/en/latest/content/tools/multiinter.html) for cohort sample-count intersections.

In a packaged implementation, Python should orchestrate `mosdepth` and `bedtools`, but it should own the simple tabular transformations. That means no required AWK, `cut`, or other shell text-processing dependency in the final tool.

## Packaged software boundary

Keep these as external subprocess calls:

- `mosdepth`
- `bedtools sort`
- `bedtools merge`
- `bedtools intersect`
- `bedtools multiinter`

Implement these in Python with streaming file readers/writers:

- Parsing `samples.tsv`.
- Extracting `PASS` intervals from `mosdepth` quantized BED output.
- Normalizing `targets.bed` to three columns.
- Counting BED interval lengths.
- Collapsing adjacent intervals by sample count.
- Collapsing adjacent intervals by population-count vector.
- Writing summary TSV files.

## Inputs

- `samples.tsv`: sample ID, population ID, and BAM/CRAM path. If populations are not needed, use a placeholder population such as `all`.
- `threshold`: minimum depth required for a site to pass in a sample.
- Optional `targets.bed`: regions to include; normalize to sorted, three-column BED before clipping cohort tables.
- Optional filters: `mapq`, excluded flags, CRAM reference FASTA, and strict/fast depth mode.

By default, omit zero-count intervals. Any site absent from the final quantized BED is assumed to have zero passing samples.

## Step 1: make per-sample pass intervals

For each sample, run `mosdepth` in quantized mode with a fail/pass split at the requested threshold:

```bash
MOSDEPTH_Q0=FAIL \
MOSDEPTH_Q1=PASS \
mosdepth \
  --threads "${threads}" \
  --fast-mode \
  --no-per-base \
  --quantize "0:${threshold}:" \
  "work/${sample}.d${threshold}" \
  "${bam_or_cram}"
```

Then keep only intervals where the sample depth is at least the threshold. In the packaged tool, implement this as a streaming Python function that reads the bgzipped `mosdepth` quantized BED and writes a three-column pass BED.

Prototype shell equivalent: read `${sample}.d${threshold}.quantized.bed.gz`, keep rows where column 4 is `PASS`, write `chrom`, `start`, and `end`, then run `bedtools sort` and `bedtools merge`.

Implementation responsibility:

- `extract_pass_intervals(quantized_bed_gz, pass_bed, pass_label="PASS")`
- `sort_and_merge_bed(pass_bed, merged_pass_bed)` using `bedtools sort` and `bedtools merge`

If `targets.bed` is provided, clip each pass BED to those targets before downstream counting:

```bash
bedtools intersect \
  -a "work/${sample}.d${threshold}.pass.bed" \
  -b "${targets_bed}" |
  bedtools sort -i - |
  bedtools merge -i - \
  > "work/${sample}.d${threshold}.pass.targets.bed"
```

## Step 2: count passing sites per sample

For each per-sample pass BED, sum interval lengths in Python. If target clipping was used in Step 1, count the target-clipped pass BED instead.

Implementation responsibility:

- `count_bed_sites(pass_bed) -> int`
- Write `results/sample_sites_at_depth_${threshold}.tsv` with columns `sample_id` and `passing_sites`.

This reports the number of genomic sites in each BAM that clear the threshold.

## Step 3: intersect pass intervals across samples

Run `bedtools multiinter` across all per-sample pass BEDs. Use the target-clipped pass BEDs if `targets.bed` was supplied. Use `-names` so the per-sample indicator columns are stable and readable.

```bash
bedtools multiinter \
  -header \
  -names ${sample_names} \
  -i ${pass_beds} \
  > "work/cohort.d${threshold}.multiinter.tsv"
```

Important output columns:

- `chrom`, `start`, `end`: sub-interval coordinates.
- `num`: number of input BEDs whose pass intervals overlap the sub-interval. This is always greater than zero in the default output.
- `list`: sample names or file IDs that pass in the sub-interval.
- One 0/1 column per sample, indicating whether that sample clears the depth threshold in the sub-interval.

Without `-empty`, `multiinter` only reports intervals overlapped by at least one input file. This is the desired behavior: intervals not present in the output are interpreted as zero passing samples.

Implementation responsibility:

- `normalize_targets_bed(targets_bed, targets_3col_bed)` to write only `chrom`, `start`, and `end` before per-sample target clipping.

## Step 4: collapse adjacent intervals with the same count

For a count-only quantized BED, keep `chrom`, `start`, `end`, and `num`, then merge adjacent intervals if they have the same count. Implement this as a streaming Python collapse over the sorted `multiinter` table.

Implementation responsibility:

- `collapse_sample_counts(multiinter_tsv, output_bed)`
- Skip the header.
- Track the current `chrom`, `start`, `end`, and count.
- Extend the current interval when the next row is adjacent and has the same count.
- Flush when the count, contig, or adjacency changes.

When `targets.bed` is supplied, this same `multiinter` table is already target-limited because the per-sample pass BEDs were clipped before `multiinter`.

This produces a `mosdepth`-like quantized file where the fourth column is the number of samples meeting the depth threshold across that interval.

If exact sample membership matters, collapse only when both `num` and the per-sample 0/1 vector are identical. That output is larger, but it preserves which samples passed rather than only how many passed.

## Step 5: collapse intervals by population counts

For multiple populations, use the per-sample 0/1 columns from `multiinter` to sum passing samples within each population. Then collapse adjacent intervals only when the full population-count vector is identical.

Assume `samples.tsv` has three tab-delimited columns:

```text
sample_id  population  bam_or_cram
```

The `sample_id` values must match the names supplied to `bedtools multiinter -names`.

CLI-level command shape:

```bash
sprite-mask collapse-populations \
  --samples samples.tsv \
  --multiinter "work/cohort.d${threshold}.multiinter.tsv" \
  --output "results/cohort.d${threshold}.population_count_quantized.bed"
```

When `targets.bed` is supplied, use the same `multiinter` table; it is already target-limited because the per-sample pass BEDs were clipped before `multiinter`.

The script should:

1. Read sample-to-population assignments from `samples.tsv`.
2. Read the `multiinter` header to identify the per-sample indicator columns.
3. For each interval, sum the 0/1 sample indicators within each population.
4. Emit `chrom`, `start`, `end`, then one count column per population.
5. Collapse adjacent intervals when `chrom`, adjacency, and all population counts match.

Concrete implementation sketch:

```python
#!/usr/bin/env python3
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--multiinter", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_population_map(path):
    sample_to_pop = {}
    populations = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "sample_id" and row[1] == "population":
                continue
            sample, population = row[0], row[1]
            sample_to_pop[sample] = population
            if population not in populations:
                populations.append(population)
    return sample_to_pop, populations


def flush(handle, current):
    if current is not None:
        chrom, start, end, counts = current
        handle.write("\t".join([chrom, str(start), str(end), *map(str, counts)]) + "\n")


def main():
    args = parse_args()
    sample_to_pop, populations = read_population_map(args.samples)
    current = None

    with open(args.multiinter) as source, open(args.output, "w") as out:
        header = source.readline().rstrip("\n").split("\t")
        sample_columns = header[5:]
        pop_index = {pop: i for i, pop in enumerate(populations)}
        column_to_pop = [sample_to_pop[sample] for sample in sample_columns]

        out.write("\t".join(["chrom", "start", "end", *populations]) + "\n")

        for line in source:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            indicators = fields[5:]
            counts = [0] * len(populations)

            for value, population in zip(indicators, column_to_pop):
                counts[pop_index[population]] += int(value)

            if current and current[0] == chrom and current[2] == start and current[3] == counts:
                current = (current[0], current[1], end, current[3])
            else:
                flush(out, current)
                current = (chrom, start, end, counts)

        flush(out, current)


if __name__ == "__main__":
    main()
```

Example output:

```text
chrom  start  end    popA  popB  popC
chr1   0      120    8     6     7
chr1   120    400    8     5     7
chr1   400    900    7     5     7
```

This is the population-aware equivalent of the count-only quantized BED: each interval is maximal with respect to an unchanged vector of passing sample counts per population.

## Step 6: summarize cohort callable space

From the quantized BED, generate a compact count distribution in Python.

Implementation responsibility:

- `summarize_sample_count_bed(sample_count_bed, output_tsv)`
- For each interval, add `end - start` to the bucket for count column 4.
- Write columns `passing_sample_count` and `sites`.

This gives the number of sites covered by exactly `N` samples at the requested depth threshold. Downstream filters can then choose, for example, regions covered by at least 80 percent of samples.

For the population-aware BED, summarize each population separately:

```bash
sprite-mask summarize-populations \
  --population-count-bed "results/cohort.d${threshold}.population_count_quantized.bed" \
  --output-prefix "results/cohort.d${threshold}.population_sites_by_count"
```

The summary should report, for each population, the number of sites covered by exactly `N` samples from that population.

## Notes and caveats

- All BED inputs to `multiinter` must be sorted in the same order.
- BED intervals are 0-based, half-open; site counts are interval lengths, `end - start`.
- Python should stream BED/TSV inputs line by line. Do not load whole genome-wide interval files into memory.
- The count-only collapse can merge adjacent regions where different samples pass, as long as the same number of samples pass. That is correct for cohort sample-count masks, but not for per-sample membership masks.
- The population-aware collapse can also merge adjacent regions where different individuals pass, as long as each population has the same passing-sample count. Preserve per-sample membership instead if the identity of samples matters.
- `mosdepth --quantize "0:${threshold}:"` treats the upper endpoint as non-inclusive, so the `PASS` bin is depth `>= threshold`.
- For `threshold=1`, the pass set is all sites with any coverage. For `threshold=0`, skip `mosdepth` and treat the full genome or target BED as passing for every sample.
- Keep the raw per-sample quantized BEDs for QC and reproducibility.
