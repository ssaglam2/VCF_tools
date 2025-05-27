# VCF Processing Pipeline: Merge, Chromosome Filter, and PASS Filter

**Author:** ssaglam2
**Date:** 27.05.2025
**Script Version:** 1.5

## Overview

This Bash script provides a flexible pipeline for processing multiple single-sample VCF (Variant Call Format) files. It automates the following key steps:

1.  **Preparation of Individual VCFs:**
    *   Identifies `.vcf.gz` files in a specified input directory.
    *   Optionally filters these files based on a user-defined substring in their filenames (e.g., to process only files like `*SAMPLE*.hard-filtered.vcf.gz`).
    *   Renames the internal sample ID within each VCF to a new ID derived from its filename, ensuring traceability.
    *   Indexes the input and reheadered VCF files.
2.  **Merging:**
    *   Merges all successfully prepared single-sample VCFs into a single multi-sample cohort VCF. This initial merge includes all variants from all chromosomes present in the input files.
3.  **Chromosome Filtering (Optional):**
    *   Filters the merged cohort VCF to retain variants only from a user-specified list of chromosomes (e.g., autosomes 1-22). This step is skipped if no chromosomes are specified.
4.  **PASS Variant Filtering:**
    *   Filters the (potentially chromosome-filtered) VCF to keep only variants that have "PASS" in their FILTER column. This is a common step to select high-quality variants.
5.  **Output:**
    *   Generates intermediate reheadered VCFs, a fully merged VCF, an optionally chromosome-filtered VCF, and a final PASS-filtered VCF, all in compressed format (`.vcf.gz`) with corresponding indexes (`.tbi`).

## Prerequisites

*   **Bash Shell:** The script is designed for a Bash environment (common on Linux, macOS, and available on Windows via WSL - Windows Subsystem for Linux).
*   **bcftools:** Version 1.9 or later is highly recommended. `bcftools` (and its dependency `htslib` which provides `bgzip` and `tabix`) must be installed and accessible in your system's `PATH`.
    *   The easiest way to install `bcftools` is often via Conda:
        ```bash
        conda create --name vcf_env -c bioconda bcftools
        conda activate vcf_env
        ```

## Configuration

Before running the script, you **must** review and potentially modify the configuration variables located at the top of the `run_vcf_processing.sh` (or your chosen script name) file.

```bash
# ====================================================================================
# --- CONFIGURATION - USER TO MODIFY THESE VALUES AS NEEDED ---
# ====================================================================================
# 1. Directory containing the input single-sample VCF.gz files.
INPUT_VCF_DIR="."

# 2. Subdirectory name to store intermediate reheadered VCF files.
REHEADERED_SUBDIR="reheadered_vcfs_intermediate"

# 3. Expected filename substring (optional).
EXPECTED_FILENAME_SUBSTRING=".hard-filtered.vcf.gz"

# 4. Suffix to strip from filename to create the sample name (override, optional).
BASENAME_STRIP_SUFFIX_OVERRIDE=""

# 5. Basename for the output merged VCF files.
OUTPUT_BASENAME="cohort_pipeline_merged"

# 6. Chromosomes to keep.
CHROMS_TO_KEEP="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
# ====================================================================================
