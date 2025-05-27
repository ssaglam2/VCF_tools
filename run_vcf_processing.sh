#!/bin/bash

# ####################################################################################
# Flexible VCF Merging, Chromosome Filtering (1-22), and PASS Filtering Pipeline
#
# Script Version: 1.5
# Author: ssaglam2
# Date: 27.05.2025
# 
#
# Description:
# This script automates:
#   1. Identifying and preparing single-sample VCF files.
#   2. Merging them into a multi-sample cohort VCF.
#   3. Filtering the merged VCF to keep only specified chromosomes (defaulting to 1-22).
#   4. Filtering the chromosome-filtered VCF for "PASS" variants.
#
# Prerequisites:
#   - Bash shell (common on Linux and macOS; WSL on Windows)
#   - bcftools (version 1.9 or later recommended) installed and accessible in the system PATH.
#     (Typically installed via Conda: `conda install -c bioconda bcftools`)
#
#
# ####################################################################################

# --- Script Behavior Options ---
set -e # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Exit status of a pipeline is that of the last command to exit with non-zero status.

# ====================================================================================
# --- CONFIGURATION - USER TO MODIFY THESE VALUES AS NEEDED ---
# ====================================================================================
# 1. Directory containing the input single-sample VCF.gz files.
#    Default is the current directory where the script is run.
INPUT_VCF_DIR="."

# 2. Subdirectory name to store intermediate reheadered VCF files.
REHEADERED_SUBDIR="reheadered_vcfs_intermediate"

# 3. Expected filename substring (optional).
#    If set, only .vcf.gz files *containing* this substring will be processed.
#    If empty (""), ALL .vcf.gz files in INPUT_VCF_DIR will be considered candidates.
#    Example for your files: ".hard-filtered.vcf.gz"
EXPECTED_FILENAME_SUBSTRING=".hard-filtered.vcf.gz"

# 4. Suffix to strip from filename to create the sample name (override, optional).
#    See script logic below for how ACTUAL_BASENAME_STRIP_SUFFIX is determined.
#    This override is for advanced cases. For your current need, leaving it empty is fine.
BASENAME_STRIP_SUFFIX_OVERRIDE=""

# 5. Basename for the output merged VCF files.
OUTPUT_BASENAME="cohort_pipeline_merged"

# 6. Chromosomes to keep.
#    Space-separated list. IMPORTANT: Must match chromosome naming in your VCFs (e.g., "1" or "chr1").
CHROMS_TO_KEEP="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
# ====================================================================================
# --- SCRIPT LOGIC --- (Generally, no need to modify below this line)
# ====================================================================================

# --- Initial Checks ---
if ! command -v bcftools &> /dev/null; then
    echo "[FATAL ERROR] bcftools command not found. Please install bcftools and ensure it is in your PATH."
    exit 1
fi
BCFTOOLS_VERSION=$(bcftools --version | head -n 1)
echo "[INFO] bcftools found: ${BCFTOOLS_VERSION}"


# Determine the actual suffix to strip for basename
ACTUAL_BASENAME_STRIP_SUFFIX=""
if [ -n "${BASENAME_STRIP_SUFFIX_OVERRIDE}" ]; then
    ACTUAL_BASENAME_STRIP_SUFFIX="${BASENAME_STRIP_SUFFIX_OVERRIDE}"
elif [ -n "${EXPECTED_FILENAME_SUBSTRING}" ]; then
    ACTUAL_BASENAME_STRIP_SUFFIX="${EXPECTED_FILENAME_SUBSTRING}"
else
    ACTUAL_BASENAME_STRIP_SUFFIX=".vcf.gz" # Default if no other suffix info given
fi

# --- Display Effective Configuration ---
echo "======================================================"
echo " VCF Merging, Chromosome (1-22) & PASS Filtering Pipeline"
echo "======================================================"
echo "Effective Configuration:"
echo "  Input VCF Directory: ${INPUT_VCF_DIR}"
echo "  Intermediate Reheadered Dir: ${REHEADERED_SUBDIR}"
if [ -n "${EXPECTED_FILENAME_SUBSTRING}" ]; then
    echo "  Processing only .vcf.gz files containing: '${EXPECTED_FILENAME_SUBSTRING}'"
else
    echo "  Processing ALL .vcf.gz files in the input directory."
fi
echo "  Suffix stripped for sample names: '${ACTUAL_BASENAME_STRIP_SUFFIX}'"
echo "  Output Merged File Basename: '${OUTPUT_BASENAME}'"
if [ -n "${CHROMS_TO_KEEP}" ]; then
    echo "  Keeping only chromosomes: ${CHROMS_TO_KEEP}"
else
    echo "  Keeping all chromosomes (CHROMS_TO_KEEP is empty)."
fi
echo "------------------------------------------------------"

# --- Main Script ---
# Navigate to the input VCF directory
cd "${INPUT_VCF_DIR}"
if [ $? -ne 0 ]; then
  echo "[ERROR] Could not change to directory '${INPUT_VCF_DIR}'. Exiting."
  exit 1
fi
ABS_INPUT_VCF_DIR=$(pwd) # Get absolute path
echo "[INFO] Working Directory: ${ABS_INPUT_VCF_DIR}"

# Create subdirectory for reheadered files
mkdir -p "${REHEADERED_SUBDIR}"
ABS_REHEADERED_SUBDIR="${ABS_INPUT_VCF_DIR}/${REHEADERED_SUBDIR}" # Absolute path

FILES_TO_MERGE_ARRAY=()
CANDIDATE_FILES_COUNT=0
PROCESSED_FILES_COUNT=0
FAILED_INDIVIDUAL_PROCESSING_COUNT=0

echo "[INFO] Stage 1: Preparing and Reheadering Individual VCFs..."
shopt -s nullglob # If no files match, loop doesn't run with literal '*'

for INPUT_VCF_GZ in *.vcf.gz; do
    CANDIDATE_FILES_COUNT=$((CANDIDATE_FILES_COUNT + 1))
    echo "  Found candidate: ${INPUT_VCF_GZ}"

    if [ -n "${EXPECTED_FILENAME_SUBSTRING}" ]; then
        if [[ "${INPUT_VCF_GZ}" != *"${EXPECTED_FILENAME_SUBSTRING}"* ]]; then
            echo "    INFO: Skipping '${INPUT_VCF_GZ}' (does not contain '${EXPECTED_FILENAME_SUBSTRING}')."
            continue
        fi
    fi

    PROCESSED_FILES_COUNT=$((PROCESSED_FILES_COUNT + 1))
    echo "  Processing target file: ${INPUT_VCF_GZ}"

    if [ ! -f "${INPUT_VCF_GZ}" ]; then # Should be redundant with nullglob but good practice
        echo "    ERROR: File '${INPUT_VCF_GZ}' is not a regular file. Skipping."
        FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
        continue
    fi

    echo "    Indexing '${INPUT_VCF_GZ}'..."
    if ! bcftools index -f -t "${INPUT_VCF_GZ}"; then
        echo "    ERROR: Failed to index ${INPUT_VCF_GZ}. Skipping."
        FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
        continue
    fi

    CURRENT_SAMPLE_NAME=$(bcftools query -l "${INPUT_VCF_GZ}")
    if [ -z "${CURRENT_SAMPLE_NAME}" ]; then
        echo "    ERROR: Could not extract sample name from ${INPUT_VCF_GZ}. Skipping."
        FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
        continue
    fi
    if [[ "${CURRENT_SAMPLE_NAME}" == *" "* ]]; then
        echo "    ERROR: ${INPUT_VCF_GZ} contains multiple sample names ('${CURRENT_SAMPLE_NAME}'). Script expects single-sample VCFs. Skipping."
        FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
        continue
    fi

    NEW_SAMPLE_NAME_RAW=$(basename "${INPUT_VCF_GZ}" "${ACTUAL_BASENAME_STRIP_SUFFIX}")
    if [ -z "${NEW_SAMPLE_NAME_RAW}" ] || [ "${NEW_SAMPLE_NAME_RAW}" == "${INPUT_VCF_GZ}" ]; then
        echo "    WARNING: Basename stripping for '${INPUT_VCF_GZ}' with suffix '${ACTUAL_BASENAME_STRIP_SUFFIX}' might not be optimal."
        NEW_SAMPLE_NAME_FALLBACK=$(basename "${INPUT_VCF_GZ}" .vcf.gz)
        if [ -z "${NEW_SAMPLE_NAME_FALLBACK}" ] || [ "${NEW_SAMPLE_NAME_FALLBACK}" == "${INPUT_VCF_GZ}" ]; then
            # Last resort: remove any final extension like .gz, then .vcf if that was it
            TMP_NAME="${INPUT_VCF_GZ%.gz}"
            NEW_SAMPLE_NAME_FALLBACK="${TMP_NAME%.vcf}"
        fi
        NEW_SAMPLE_NAME_RAW="${NEW_SAMPLE_NAME_FALLBACK}"
    fi
    NEW_SAMPLE_NAME_SANITIZED=$(echo "${NEW_SAMPLE_NAME_RAW}" | sed 's/[^a-zA-Z0-9_-]/_/g' | sed 's/^_//;s/_$//') # Sanitize and remove leading/trailing underscores

    echo "    Current internal sample: '${CURRENT_SAMPLE_NAME}' -> New sanitized sample name: '${NEW_SAMPLE_NAME_SANITIZED}'"

    REHEADERED_VCF_PATH="${ABS_REHEADERED_SUBDIR}/${NEW_SAMPLE_NAME_SANITIZED}.renamed.vcf.gz"
    RENAME_MAP_FILE_PATH="${ABS_REHEADERED_SUBDIR}/tmp_rename_map.${NEW_SAMPLE_NAME_SANITIZED}.$$" # Unique temp file with PID

    echo -e "${CURRENT_SAMPLE_NAME}\t${NEW_SAMPLE_NAME_SANITIZED}" > "${RENAME_MAP_FILE_PATH}"

    echo "    Reheadering to: ${REHEADERED_VCF_PATH}"
    if bcftools reheader -s "${RENAME_MAP_FILE_PATH}" "${INPUT_VCF_GZ}" -o "${REHEADERED_VCF_PATH}"; then
        echo "    Indexing reheadered file: ${REHEADERED_VCF_PATH}"
        if ! bcftools index -f -t "${REHEADERED_VCF_PATH}"; then
             echo "    ERROR: Failed to index reheadered file ${REHEADERED_VCF_PATH}. Skipping this file for merge."
             FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
        else
            FILES_TO_MERGE_ARRAY+=("${REHEADERED_VCF_PATH}")
        fi
    else
        echo "    ERROR: Failed to reheader ${INPUT_VCF_GZ}. Skipping."
        FAILED_INDIVIDUAL_PROCESSING_COUNT=$((FAILED_INDIVIDUAL_PROCESSING_COUNT + 1))
    fi
    rm -f "${RENAME_MAP_FILE_PATH}" # Clean up temp map file
done
shopt -u nullglob

echo "[INFO] Finished preparing individual VCFs."
echo "  Total candidates found (*.vcf.gz): ${CANDIDATE_FILES_COUNT}"
echo "  Files matching substring filter (if any) and targeted for processing: ${PROCESSED_FILES_COUNT}"
echo "  Files successfully prepared for merging: ${#FILES_TO_MERGE_ARRAY[@]}"
echo "  Files that failed during individual processing: ${FAILED_INDIVIDUAL_PROCESSING_COUNT}"
echo "------------------------------------------------------"


# --- Merging VCFs ---
MERGED_ALL_CHROMS_VCF="${OUTPUT_BASENAME}_all_chroms.vcf.gz" # This will be in INPUT_VCF_DIR (current dir)
if [ ${#FILES_TO_MERGE_ARRAY[@]} -eq 0 ]; then
    echo "[ERROR] No VCF files were successfully prepared for merging. Exiting."
    exit 1
fi

echo "[INFO] Stage 2: Merging ${#FILES_TO_MERGE_ARRAY[@]} VCF files into ${MERGED_ALL_CHROMS_VCF}..."
if bcftools merge "${FILES_TO_MERGE_ARRAY[@]}" -m none -Oz -o "${MERGED_ALL_CHROMS_VCF}"; then
    echo "  Merging successful."
    echo "  Indexing merged file: ${MERGED_ALL_CHROMS_VCF}"
    if ! bcftools index -f -t "${MERGED_ALL_CHROMS_VCF}"; then
        echo "  [WARNING] Failed to index the merged file ${MERGED_ALL_CHROMS_VCF}. Further steps might be affected."
    else
        echo "  Merged file indexed successfully."
    fi
else
    echo "[ERROR] bcftools merge command failed. Exiting."
    exit 1
fi
echo "------------------------------------------------------"


# This will be the input for the PASS filter stage
CURRENT_MERGED_VCF_TO_FILTER="${MERGED_ALL_CHROMS_VCF}"
MERGED_CHROMS_FILTERED_VCF="" # Initialize to empty

# --- STAGE 2.5: Chromosome Filtering ---
if [ -n "${CHROMS_TO_KEEP}" ]; then
    # Create a more concise filename part for chromosome filtered file
    CHROMS_LABEL_PART_RAW=$(echo ${CHROMS_TO_KEEP} | tr ' ' '_')
    # Further sanitize and shorten if necessary for filename
    CHROMS_LABEL_PART=$(echo "${CHROMS_LABEL_PART_RAW}" | sed 's/[^a-zA-Z0-9_]/_/g' | cut -c1-30)
    if [ "${CHROMS_LABEL_PART_RAW}" == "1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22" ]; then
      CHROMS_LABEL_PART="autosomes_1_22" # Specific label for common case
    elif [ ${#CHROMS_LABEL_PART_RAW} -gt 30 ]; then
      CHROMS_LABEL_PART="selected_chroms"
    fi
    MERGED_CHROMS_FILTERED_VCF="${OUTPUT_BASENAME}_${CHROMS_LABEL_PART}.vcf.gz"
    
    echo "[INFO] Stage 2.5: Filtering for specified chromosomes into ${MERGED_CHROMS_FILTERED_VCF}..."

    if [ ! -f "${MERGED_ALL_CHROMS_VCF}" ]; then # Should exist from previous stage
        echo "[ERROR] File ${MERGED_ALL_CHROMS_VCF} not found for chromosome filtering. Exiting."
        exit 1
    fi
    # Ensure index for merged_all_chroms.vcf.gz exists
    if [ ! -f "${MERGED_ALL_CHROMS_VCF}.tbi" ] && [ ! -f "${MERGED_ALL_CHROMS_VCF}.csi" ]; then
        echo "[WARNING] Index for ${MERGED_ALL_CHROMS_VCF} not found. Attempting to index..."
        if ! bcftools index -f -t "${MERGED_ALL_CHROMS_VCF}"; then
            echo "[ERROR] Failed to re-index ${MERGED_ALL_CHROMS_VCF}. Cannot filter chromosomes. Exiting."
            exit 1
        fi
    fi

    # Convert space-separated CHROMS_TO_KEEP to comma-separated for bcftools regions
    REGIONS_TO_KEEP=$(echo "${CHROMS_TO_KEEP}" | tr ' ' ',')
    echo "  Filtering for regions: ${REGIONS_TO_KEEP}"

    if bcftools view --regions "${REGIONS_TO_KEEP}" "${MERGED_ALL_CHROMS_VCF}" -Oz -o "${MERGED_CHROMS_FILTERED_VCF}"; then
        echo "  Chromosome filtering successful."
        echo "  Indexing chromosome-filtered file: ${MERGED_CHROMS_FILTERED_VCF}"
        if ! bcftools index -f -t "${MERGED_CHROMS_FILTERED_VCF}"; then
            echo "  [WARNING] Failed to index ${MERGED_CHROMS_FILTERED_VCF}."
        else
            echo "  Chromosome-filtered file indexed successfully."
        fi
        CURRENT_MERGED_VCF_TO_FILTER="${MERGED_CHROMS_FILTERED_VCF}" # Update for next stage
    else
        echo "[ERROR] Chromosome filtering failed. PASS filtering will proceed on the all-chromosome merged file."
        # CURRENT_MERGED_VCF_TO_FILTER remains MERGED_ALL_CHROMS_VCF
        MERGED_CHROMS_FILTERED_VCF="" # Unset as it failed or wasn't created
    fi
    echo "------------------------------------------------------"
else
    echo "[INFO] Stage 2.5: Chromosome filtering skipped as CHROMS_TO_KEEP is not set."
    echo "------------------------------------------------------"
fi


# --- Filtering Merged VCF for PASS variants ---
MERGED_PASS_VCF="${OUTPUT_BASENAME}_final_PASS.vcf.gz" # More descriptive final name
echo "[INFO] Stage 3: Filtering ${CURRENT_MERGED_VCF_TO_FILTER} for PASS variants -> ${MERGED_PASS_VCF}..."

if [ ! -f "${CURRENT_MERGED_VCF_TO_FILTER}" ]; then
    echo "[ERROR] Input VCF file for PASS filtering ('${CURRENT_MERGED_VCF_TO_FILTER}') not found. Exiting."
    exit 1
fi
# Ensure index for the VCF to be PASS-filtered
if [ ! -f "${CURRENT_MERGED_VCF_TO_FILTER}.tbi" ] && [ ! -f "${CURRENT_MERGED_VCF_TO_FILTER}.csi" ]; then
    echo "[WARNING] Index for ${CURRENT_MERGED_VCF_TO_FILTER} not found. Attempting to index..."
    if ! bcftools index -f -t "${CURRENT_MERGED_VCF_TO_FILTER}"; then
        echo "[ERROR] Failed to index ${CURRENT_MERGED_VCF_TO_FILTER}. Cannot PASS filter. Exiting."
        exit 1
    fi
fi

if bcftools view --include 'FILTER="PASS"' "${CURRENT_MERGED_VCF_TO_FILTER}" -Oz -o "${MERGED_PASS_VCF}"; then
    echo "  PASS filtering successful."
    echo "  Indexing PASS-filtered file: ${MERGED_PASS_VCF}"
    if ! bcftools index -f -t "${MERGED_PASS_VCF}"; then
        echo "  [WARNING] Failed to index the PASS-filtered file ${MERGED_PASS_VCF}."
    else
        echo "  PASS-filtered file indexed successfully."
    fi
else
    echo "[ERROR] bcftools view for PASS filtering failed. The input file (${CURRENT_MERGED_VCF_TO_FILTER}) is still available."
fi
echo "------------------------------------------------------"

echo "[INFO] Pipeline Finished!"
echo "  Input VCF Directory processed: ${ABS_INPUT_VCF_DIR}"
echo "  Intermediate reheadered files are in: ${ABS_REHEADERED_SUBDIR}/"
echo "  Initial merged VCF (all chroms): ${ABS_INPUT_VCF_DIR}/${MERGED_ALL_CHROMS_VCF}"
if [ -n "${MERGED_CHROMS_FILTERED_VCF}" ] && [ -f "${MERGED_CHROMS_FILTERED_VCF}" ]; then # Check if var is set and file exists
    echo "  Merged VCF (filtered for specific chroms): ${ABS_INPUT_VCF_DIR}/${MERGED_CHROMS_FILTERED_VCF}"
fi
if [ -f "${MERGED_PASS_VCF}" ]; then
    echo "  Final VCF (PASS variants from '${CURRENT_MERGED_VCF_TO_FILTER}'): ${ABS_INPUT_VCF_DIR}/${MERGED_PASS_VCF}"
else
    echo "  PASS-filtered VCF was not successfully created."
fi

echo ""
echo "Summary Counts:"
if [ -f "${MERGED_ALL_CHROMS_VCF}" ]; then
    RAW_COUNT=$(bcftools view -H "${MERGED_ALL_CHROMS_VCF}" | grep -vc '^#' || echo "N/A")
    echo "  Variants in initial merged file (${MERGED_ALL_CHROMS_VCF}): ${RAW_COUNT}"
fi
if [ -n "${MERGED_CHROMS_FILTERED_VCF}" ] && [ -f "${MERGED_CHROMS_FILTERED_VCF}" ]; then
    CHROMS_FILTERED_COUNT=$(bcftools view -H "${MERGED_CHROMS_FILTERED_VCF}" | grep -vc '^#' || echo "N/A")
    echo "  Variants after chromosome filtering (${MERGED_CHROMS_FILTERED_VCF}): ${CHROMS_FILTERED_COUNT}"
fi
if [ -f "${MERGED_PASS_VCF}" ]; then
    PASS_COUNT=$(bcftools view -H "${MERGED_PASS_VCF}" | grep -vc '^#' || echo "N/A")
    echo "  Variants in final PASS-filtered file (${MERGED_PASS_VCF}): ${PASS_COUNT}"
fi
echo "======================================================"
