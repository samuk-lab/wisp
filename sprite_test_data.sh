#!/usr/bin/env bash
set -euo pipefail

SCRIPT_VERSION="2026-05-25-sprite-test-data-env-v4-single-samtools-view"

# download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#
# macOS/Linux-compatible script to:
#   1. Create a conda environment using mamba/micromamba/libmamba only.
#   2. Download a small high-coverage 1000 Genomes two-population dataset.
#
# This version DOES NOT use `conda activate`. It creates the environment at an
# explicit prefix and runs tools with `conda run -p` or `micromamba run -p`.
# This avoids EnvironmentNameNotFound errors caused by mismatched conda/mamba
# environment registries.
#
# Dataset:
#   - 10 regional BAMs: 5 GBR + 5 YRI
#   - one indexed multi-sample VCF containing all 10 samples
#   - samples.tsv metadata for sprite-mask
#   - targets.bed for the selected region
#
# Source:
#   - Alignments: 1000 Genomes / NYGC ~30x GRCh38 CRAMs
#   - VCF: 1000 Genomes high-coverage GRCh38 chr20 multi-sample VCF
#
# By default, regional BAMs are downsampled from ~30x to roughly ~20x:
#   DOWNSAMPLE_FRAC=0.67
#
# Set DOWNSAMPLE_FRAC=1 to keep full ~30x regional BAMs.
#
# Requirements:
#   - macOS or Linux
#   - bash
#   - conda+mamba, micromamba, or conda with --solver libmamba
#
# Usage:
#   bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#
# Optional:
#   ENV_NAME=sprite-test-data bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#   OUTDIR=tests/test_data/1000g_10sample_highcov_subset bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#   REGION=chr20:10000000-10100000 bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#   THREADS=4 bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#   DOWNSAMPLE_FRAC=1 bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#   FORCE=1 bash download_1000g_10sample_highcov_subset_sprite_test_data_v4.sh
#
# Notes:
#   - REGION must be on chr20 unless you also update VCF_URL.
#   - REGION uses samtools/bcftools coordinates: 1-based inclusive.
#   - targets.bed is written as BED: 0-based half-open.
#   - Source CRAMs and VCF are GRCh38 / chr-prefixed.
#   - Local alignments are written as BAMs so downstream tests do not need
#     CRAM reference handling.

ENV_NAME="${ENV_NAME:-sprite-test-data}"
OUTDIR="${OUTDIR:-tests/test_data/1000g_10sample_highcov_subset}"
REGION="${REGION:-chr20:10000000-10100000}"
THREADS="${THREADS:-2}"
FORCE="${FORCE:-0}"

# Set to 1 for full ~30x. Default 0.67 gives roughly 20x from 30x source data.
DOWNSAMPLE_FRAC="${DOWNSAMPLE_FRAC:-0.67}"
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-42}"

# 1000 Genomes high-coverage 3202-sample phased SNV/INDEL/SV VCF, GRCh38.
VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Format: sample_id population cram_url
# CRAMs are from the 1000G 2504 high-coverage sequence index.
SAMPLES=(
  "HG00096 GBR https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram"
  "HG00097 GBR https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240115/HG00097.final.cram"
  "HG00099 GBR https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240116/HG00099.final.cram"
  "HG00100 GBR https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240117/HG00100.final.cram"
  "HG00101 GBR https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240118/HG00101.final.cram"
  "NA18486 YRI https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239335/NA18486.final.cram"
  "NA18489 YRI https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239336/NA18489.final.cram"
  "NA18498 YRI https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239337/NA18498.final.cram"
  "NA18499 YRI https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239338/NA18499.final.cram"
  "NA18501 YRI https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239339/NA18501.final.cram"
)

die() {
    echo "ERROR: $*" >&2
    exit 1
}

info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_os() {
    case "$(uname -s)" in
        Darwin|Linux)
            ;;
        *)
            die "This script supports macOS and Linux only."
            ;;
    esac
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Global variables set by create_env_and_runner:
ENV_PREFIX=""
RUNNER_KIND=""
CREATE_TOOL=""

create_env_and_runner() {
    info "Script version: ${SCRIPT_VERSION}"
    info "Environment name: ${ENV_NAME}"

    # Preferred: conda + mamba. Create under conda base envs so the prefix is predictable.
    if command_exists conda && command_exists mamba; then
        CONDA_BASE="$(conda info --base)"
        ENV_PREFIX="${CONDA_BASE}/envs/${ENV_NAME}"
        CREATE_TOOL="mamba"
        RUNNER_KIND="conda"

        if [ ! -d "${ENV_PREFIX}" ]; then
            info "Creating environment with mamba at: ${ENV_PREFIX}"
            mamba create -y \
                -p "${ENV_PREFIX}" \
                -c conda-forge \
                -c bioconda \
                python=3.11 \
                curl \
                samtools \
                bcftools \
                htslib \
                mosdepth \
                bedtools
        else
            info "Environment prefix already exists: ${ENV_PREFIX}"
        fi
        return
    fi

    # Alternative: micromamba only. Use an explicit prefix under the current project.
    if command_exists micromamba; then
        ENV_PREFIX="${PWD}/.conda_envs/${ENV_NAME}"
        CREATE_TOOL="micromamba"
        RUNNER_KIND="micromamba"

        if [ ! -d "${ENV_PREFIX}" ]; then
            info "Creating environment with micromamba at: ${ENV_PREFIX}"
            micromamba create -y \
                -p "${ENV_PREFIX}" \
                -c conda-forge \
                -c bioconda \
                python=3.11 \
                curl \
                samtools \
                bcftools \
                htslib \
                mosdepth \
                bedtools
        else
            info "Environment prefix already exists: ${ENV_PREFIX}"
        fi
        return
    fi

    # Fallback: conda with libmamba solver explicitly.
    if command_exists conda; then
        if ! conda create --help 2>/dev/null | grep -q -- "--solver"; then
            die "conda is installed, but it does not expose --solver. Install mamba/micromamba or a conda version with the libmamba solver."
        fi

        CONDA_BASE="$(conda info --base)"
        ENV_PREFIX="${CONDA_BASE}/envs/${ENV_NAME}"
        CREATE_TOOL="conda --solver libmamba"
        RUNNER_KIND="conda"

        if [ ! -d "${ENV_PREFIX}" ]; then
            info "Creating environment with conda --solver libmamba at: ${ENV_PREFIX}"
            conda create --solver libmamba -y \
                -p "${ENV_PREFIX}" \
                -c conda-forge \
                -c bioconda \
                python=3.11 \
                curl \
                samtools \
                bcftools \
                htslib \
                mosdepth \
                bedtools
        else
            info "Environment prefix already exists: ${ENV_PREFIX}"
        fi
        return
    fi

    die "mamba, micromamba, or conda with --solver libmamba is required."
}

run_in_env() {
    if [ "${RUNNER_KIND}" = "conda" ]; then
        conda run --no-capture-output -p "${ENV_PREFIX}" "$@"
    elif [ "${RUNNER_KIND}" = "micromamba" ]; then
        micromamba run -p "${ENV_PREFIX}" "$@"
    else
        die "internal error: RUNNER_KIND is not set"
    fi
}

validate_tools() {
    run_in_env curl --version >/dev/null
    run_in_env samtools --version >/dev/null
    run_in_env bcftools --version >/dev/null
    run_in_env mosdepth --version >/dev/null
    run_in_env bedtools --version >/dev/null

    info "Environment prefix: ${ENV_PREFIX}"
    info "Environment created with: ${CREATE_TOOL}"
    info "samtools: $(run_in_env samtools --version | head -n 1)"
    info "bcftools: $(run_in_env bcftools --version | head -n 1)"
}

validate_region_and_write_targets() {
    chrom="${REGION%%:*}"
    range="${REGION#*:}"
    start="${range%-*}"
    end="${range#*-}"

    if [ "${chrom}" != "chr20" ]; then
        die "default VCF_URL is for chr20. Use REGION=chr20:start-end or edit VCF_URL."
    fi

    case "${start}" in
        ''|*[!0-9]*) die "REGION start must be numeric, e.g. chr20:10000000-10100000" ;;
    esac

    case "${end}" in
        ''|*[!0-9]*) die "REGION end must be numeric, e.g. chr20:10000000-10100000" ;;
    esac

    if [ "${start}" -ge "${end}" ]; then
        die "REGION start must be less than REGION end."
    fi

    mkdir -p "${OUTDIR}"
    bed_start=$((start - 1))
    printf "%s\t%s\t%s\n" "${chrom}" "${bed_start}" "${end}" > "${OUTDIR}/targets.bed"
}

write_metadata_files() {
    mkdir -p "${OUTDIR}/bams"

    {
        printf "sample_id\tpopulation\talignment\n"
        for entry in "${SAMPLES[@]}"; do
            read -r sample pop cram_url <<< "${entry}"
            printf "%s\t%s\t%s\n" "${sample}" "${pop}" "${OUTDIR}/bams/${sample}.bam"
        done
    } > "${OUTDIR}/samples.tsv"

    {
        for entry in "${SAMPLES[@]}"; do
            read -r sample pop cram_url <<< "${entry}"
            printf "%s\n" "${sample}"
        done
    } > "${OUTDIR}/samples.list"

    {
        printf "sample_id\tpopulation\tsource_cram\n"
        for entry in "${SAMPLES[@]}"; do
            read -r sample pop cram_url <<< "${entry}"
            printf "%s\t%s\t%s\n" "${sample}" "${pop}" "${cram_url}"
        done
    } > "${OUTDIR}/sample_populations_and_sources.tsv"
}

configure_cram_reference_cache() {
    mkdir -p "${OUTDIR}/ref_cache"

    # htslib can retrieve CRAM reference slices by MD5 from ENA if they are not
    # already present locally. The local cache prevents repeated downloads.
    export REF_CACHE="${OUTDIR}/ref_cache/%2s/%2s/%s"
    export REF_PATH="${REF_CACHE}:https://www.ebi.ac.uk/ena/cram/md5/%s"
}

samtools_downsample_arg() {
    frac="$1"

    case "${frac}" in
        1|1.0|1.00)
            printf ""
            ;;
        0.*)
            printf "%s.%s" "${DOWNSAMPLE_SEED}" "${frac#0.}"
            ;;
        .*)
            printf "%s.%s" "${DOWNSAMPLE_SEED}" "${frac#.}"
            ;;
        *)
            die "DOWNSAMPLE_FRAC must be between 0 and 1, e.g. 0.67 or 1"
            ;;
    esac
}

download_bam_region() {
    sample="$1"
    pop="$2"
    cram_url="$3"
    crai_url="${cram_url}.crai"
    out_bam="${OUTDIR}/bams/${sample}.bam"

    if [ "${FORCE}" != "1" ] && [ -s "${out_bam}" ] && [ -s "${out_bam}.bai" ]; then
        info "BAM already exists, skipping: ${out_bam}"
        return
    fi

    info "Downloading regional BAM for ${sample} (${pop}): ${REGION}"
    info "  source: ${cram_url}"
    info "  downsample fraction: ${DOWNSAMPLE_FRAC}"

    rm -f "${out_bam}" "${out_bam}.bai"

    downsample_arg="$(samtools_downsample_arg "${DOWNSAMPLE_FRAC}")"

    if [ -z "${downsample_arg}" ]; then
        run_in_env samtools view \
            -@ "${THREADS}" \
            -b \
            -o "${out_bam}" \
            "${cram_url}##idx##${crai_url}" \
            "${REGION}"
    else
        # Do region extraction and downsampling in one samtools invocation.
        # This avoids piping between two `conda run` processes, which can fail
        # with "[main_samview] fail to read the header from '-'."
        run_in_env samtools view \
            -@ "${THREADS}" \
            -b \
            -s "${downsample_arg}" \
            -o "${out_bam}" \
            "${cram_url}##idx##${crai_url}" \
            "${REGION}"
    fi

    run_in_env samtools quickcheck -v "${out_bam}"
    run_in_env samtools index -@ "${THREADS}" "${out_bam}"
}

download_multisample_vcf_region() {
    safe_region="${REGION/:/_}"
    safe_region="${safe_region//:/_}"
    out_vcf="${OUTDIR}/1000g_10samples_highcov.${safe_region}.vcf.gz"

    if [ "${FORCE}" != "1" ] && [ -s "${out_vcf}" ] && [ -s "${out_vcf}.tbi" ]; then
        info "VCF already exists, skipping: ${out_vcf}"
        ln -sf "$(basename "${out_vcf}")" "${OUTDIR}/1000g_10samples_highcov.vcf.gz"
        ln -sf "$(basename "${out_vcf}.tbi")" "${OUTDIR}/1000g_10samples_highcov.vcf.gz.tbi"
        return
    fi

    info "Downloading multi-sample high-coverage VCF subset: ${REGION}"
    info "  source: ${VCF_URL}"

    rm -f "${out_vcf}" "${out_vcf}.tbi"

    run_in_env bcftools view \
        -r "${REGION}" \
        -S "${OUTDIR}/samples.list" \
        -Oz \
        -o "${out_vcf}" \
        "${VCF_URL}"

    run_in_env bcftools index -t "${out_vcf}"

    ln -sf "$(basename "${out_vcf}")" "${OUTDIR}/1000g_10samples_highcov.vcf.gz"
    ln -sf "$(basename "${out_vcf}.tbi")" "${OUTDIR}/1000g_10samples_highcov.vcf.gz.tbi"
}

write_dataset_readme() {
    cat > "${OUTDIR}/README.md" <<EOF
# 1000 Genomes 10-sample high-coverage test subset

This directory contains a small, regional 1000 Genomes high-coverage test dataset.

## Samples

- 5 GBR samples: HG00096, HG00097, HG00099, HG00100, HG00101
- 5 YRI samples: NA18486, NA18489, NA18498, NA18499, NA18501

## Region

\`${REGION}\`

Coordinates are GRCh38 / chr-prefixed 1000 Genomes high-coverage coordinates.

## Coverage

The source CRAMs are NYGC high-coverage, approximately 30x whole-genome data.
This script uses DOWNSAMPLE_FRAC=${DOWNSAMPLE_FRAC}, so the regional BAMs are
approximately 30x * ${DOWNSAMPLE_FRAC}. Set DOWNSAMPLE_FRAC=1 to keep full depth.

## Files

- \`samples.tsv\`: sprite-mask sample metadata
- \`sample_populations_and_sources.tsv\`: sample/population/source metadata
- \`samples.list\`: sample list used for VCF subsetting
- \`targets.bed\`: BED interval for the selected region
- \`bams/*.bam\`: regional BAMs for each sample
- \`1000g_10samples_highcov.vcf.gz\`: multi-sample VCF for all 10 samples in the selected region

## Note

The source alignments are CRAMs, but the local test files are BAMs so downstream
tests do not need to handle CRAM reference lookup.
EOF
}

main() {
    check_os
    create_env_and_runner
    validate_tools

    mkdir -p "${OUTDIR}/bams"

    validate_region_and_write_targets
    write_metadata_files
    configure_cram_reference_cache

    for entry in "${SAMPLES[@]}"; do
        read -r sample pop cram_url <<< "${entry}"
        download_bam_region "${sample}" "${pop}" "${cram_url}"
    done

    download_multisample_vcf_region
    write_dataset_readme

    echo
    info "Done."
    echo "Environment name:   ${ENV_NAME}"
    echo "Environment prefix: ${ENV_PREFIX}"
    echo "Output dir:         ${OUTDIR}"
    echo "Samples TSV:        ${OUTDIR}/samples.tsv"
    echo "Targets BED:        ${OUTDIR}/targets.bed"
    echo "Multi VCF:          ${OUTDIR}/1000g_10samples_highcov.vcf.gz"
}

main "$@"
