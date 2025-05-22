#!/bin/bash

# ####################################################################################
# Flexible VCF Merging and PASS Filtering Pipeline
#
# Script Version: 1.2
# Author: [Your Name/GitHub Username Here]
# Date: [Date you finalize this for GitHub]
#
# Description:
# This script automates the merging of multiple single-sample VCF files and
# subsequent filtering for "PASS" variants. It offers flexibility in how
# input files are identified and how sample names are derived for the merged output.
#
# Key Features:
#   - Identifies single-sample VCF files (ending in .vcf.gz).
#   - Optionally filters these candidates based on a user-defined filename substring.
#   - Renames internal VCF sample names based on (a portion of) the input filename.
#   - Merges processed VCFs into a multi-sample cohort VCF.
#   - Filters the merged cohort VCF to retain only variants with "PASS" in the FILTER column.
#   - Creates necessary output directories and indexes VCF files at each step.
#
# Prerequisites:
#   - Bash shell (common on Linux and macOS; WSL on Windows)
#   - bcftools (version 1.9 or later recommended) installed and accessible in the system PATH.
#     (Typically installed via Conda: `conda install -c bioconda bcftools`)
#
# Suggested README.md content for GitHub:
#   - Brief overview of the script.
#   - Detailed "Prerequisites" section.
#   - "Configuration" section explaining each variable at the top of the script.
#     Provide examples for `INPUT_VCF_DIR`, `EXPECTED_FILENAME_SUBSTRING`,
#     `BASENAME_STRIP_SUFFIX_OVERRIDE`, and `OUTPUT_BASENAME`.
#   - "Usage" instructions:
#     1. Clone the repository or download the script.
#     2. Modify the CONFIGURATION variables in the script to match your setup.
#     3. Make the script executable: `chmod +x run_vcf_processing.sh`
#     4. Ensure your Conda environment with bcftools is activated.
#     5. Run the script: `./run_vcf_processing.sh`
#   - "Output Files" description.
#   - "Troubleshooting" tips (e.g., check bcftools installation, file paths, permissions).
#   - "License" (e.g., MIT, GPL).
#   - "Citation" (if applicable, or how to acknowledge the script).
#
# ####################################################################################

# --- Script Behavior Options ---
set -e # Exit immediately if a command exits with a non-zero status.
# set -u # Treat unset variables as an error when substituting (can be too strict for some).
set -o pipefail # Exit status of a pipeline is that of the last command to exit with non-zero status.


# ====================================================================================
# --- CONFIGURATION - USER TO MODIFY THESE VALUES AS NEEDED ---
# ====================================================================================

# 1. Directory containing the input single-sample VCF.gz files.
#    Default is the current directory where the script is run.
#    Example: "/path/to/my/vcfs"
INPUT_VCF_DIR="."

# 2. Subdirectory name to store intermediate reheadered VCF files.
#    This will be created within the INPUT_VCF_DIR if it doesn't exist.
REHEADERED_SUBDIR="reheadered_vcfs_intermediate"

# 3. Expected filename substring (optional).
#    If set, only .vcf.gz files *containing* this substring will be processed.
#    If you want to process ALL .vcf.gz files in INPUT_VCF_DIR, set this to "" (empty string).
#    Example for files like "SAMPLE_ABC.hard-filtered.vcf.gz": ".hard-filtered.vcf.gz"
#    Example for files like "20122BioDRI.custom.vcf.gz": ".BioDRI.custom.vcf.gz"
EXPECTED_FILENAME_SUBSTRING=".hard-filtered.vcf.gz" # For your CRB/BioDRI files

# 4. Suffix to strip from filename to create the sample name (override, optional).
#    By default, if EXPECTED_FILENAME_SUBSTRING is set, that substring is used as the
#    suffix to strip by `basename`. If EXPECTED_FILENAME_SUBSTRING is empty, ".vcf.gz" is stripped.
#    If this default logic doesn't produce your desired sample ID prefix, you can explicitly
#    set BASENAME_STRIP_SUFFIX_OVERRIDE to the exact string you want `basename` to remove.
#    Example: If files are "SAMPLE_ABC.processed.vcf.gz" and you want "SAMPLE_ABC" as ID,
#             and EXPECTED_FILENAME_SUBSTRING=".processed.vcf.gz", the default is fine.
#             But if files are "SAMPLE_ABC_raw.vcf.gz" and you only set
#             EXPECTED_FILENAME_SUBSTRING="", BASENAME_STRIP_SUFFIX_OVERRIDE could be "_raw.vcf.gz".
#    For your case (stripping ".hard-filtered.vcf.gz"), the default logic works if
#    EXPECTED_FILENAME_SUBSTRING is set as above, so this can be left empty.
BASENAME_STRIP_SUFFIX_OVERRIDE=""

# 5. Basename for the output merged VCF files.
#    Example: "my_cohort" will result in "my_cohort.vcf.gz" and "my_cohort_PASS.vcf.gz".
OUTPUT_BASENAME="cohort_pipeline_merged"

# ====================================================================================
# --- SCRIPT LOGIC --- (Generally, no need to modify below this line)
# ====================================================================================

# --- Initial Checks ---
if ! command -v bcftools &> /dev/null; then
    echo "[FATAL ERROR] bcftools command not found. Please install bcftools and ensure it is in your PATH."
    exit 1
fi
echo "[INFO] bcftools found: $(bcftools --version | head -n 1)"


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
echo " VCF Merging and PASS Filtering Pipeline"
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

    # Secondary check for specific filename pattern, if defined
    if [ -n "${EXPECTED_FILENAME_SUBSTRING}" ]; then
        if [[ "${INPUT_VCF_GZ}" != *"${EXPECTED_FILENAME_SUBSTRING}"* ]]; then
            echo "    INFO: Skipping '${INPUT_VCF_GZ}' as it does not contain expected substring '${EXPECTED_FILENAME_SUBSTRING}'."
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
        echo "    WARNING: Basename stripping might not have worked as expected for '${INPUT_VCF_GZ}' with suffix '${ACTUAL_BASENAME_STRIP_SUFFIX}'."
        echo "             Using full filename without just '.vcf.gz' as fallback sample name."
        NEW_SAMPLE_NAME_RAW=$(basename "${INPUT_VCF_GZ}" .vcf.gz)
        if [ -z "${NEW_SAMPLE_NAME_RAW}" ] || [ "${NEW_SAMPLE_NAME_RAW}" == "${INPUT_VCF_GZ}" ]; then
            NEW_SAMPLE_NAME_RAW="${INPUT_VCF_GZ%.*}" # Last resort: remove any final extension
        fi
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
             # rm -f "${REHEADERED_VCF_PATH}" # Optional: clean up failed reheadered file
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
MERGED_RAW_VCF="${OUTPUT_BASENAME}.vcf.gz" # This will be in INPUT_VCF_DIR (current dir)
if [ ${#FILES_TO_MERGE_ARRAY[@]} -eq 0 ]; then
    echo "[ERROR] No VCF files were successfully prepared for merging. Exiting."
    exit 1
fi

echo "[INFO] Stage 2: Merging ${#FILES_TO_MERGE_ARRAY[@]} VCF files into ${MERGED_RAW_VCF}..."
if bcftools merge "${FILES_TO_MERGE_ARRAY[@]}" -m none -Oz -o "${MERGED_RAW_VCF}"; then
    echo "  Merging successful."
    echo "  Indexing merged file: ${MERGED_RAW_VCF}"
    if ! bcftools index -f -t "${MERGED_RAW_VCF}"; then
        echo "  [WARNING] Failed to index the merged file ${MERGED_RAW_VCF}. Further steps might be affected."
    else
        echo "  Merged file indexed successfully."
    fi
else
    echo "[ERROR] bcftools merge command failed. Exiting."
    exit 1
fi
echo "------------------------------------------------------"


# --- Filtering Merged VCF for PASS variants ---
MERGED_PASS_VCF="${OUTPUT_BASENAME}_PASS.vcf.gz" # Also in current dir
echo "[INFO] Stage 3: Filtering ${MERGED_RAW_VCF} for PASS variants -> ${MERGED_PASS_VCF}..."

if [ ! -f "${MERGED_RAW_VCF}" ]; then
    echo "[ERROR] Merged VCF file ${MERGED_RAW_VCF} not found. Cannot filter. Exiting."
    exit 1
fi
# Ensure index for merged raw VCF exists before filtering
if [ ! -f "${MERGED_RAW_VCF}.tbi" ] && [ ! -f "${MERGED_RAW_VCF}.csi" ]; then
    echo "[WARNING] Index for ${MERGED_RAW_VCF} not found. Attempting to index now..."
    if ! bcftools index -f -t "${MERGED_RAW_VCF}"; then
        echo "[ERROR] Failed to re-index ${MERGED_RAW_VCF}. Cannot filter. Exiting."
        exit 1
    fi
fi

if bcftools view --include 'FILTER="PASS"' "${MERGED_RAW_VCF}" -Oz -o "${MERGED_PASS_VCF}"; then
    echo "  PASS filtering successful."
    echo "  Indexing PASS-filtered file: ${MERGED_PASS_VCF}"
    if ! bcftools index -f -t "${MERGED_PASS_VCF}"; then
        echo "  [WARNING] Failed to index the PASS-filtered file ${MERGED_PASS_VCF}."
    else
        echo "  PASS-filtered file indexed successfully."
    fi
else
    echo "[ERROR] bcftools view for PASS filtering failed. The raw merged file is still available."
fi
echo "------------------------------------------------------"

echo "[INFO] Pipeline Finished!"
echo "  Input VCF Directory processed: ${ABS_INPUT_VCF_DIR}"
echo "  Intermediate reheadered files are in: ${ABS_REHEADERED_SUBDIR}/"
echo "  Final merged (unfiltered) VCF: ${ABS_INPUT_VCF_DIR}/${MERGED_RAW_VCF}"
if [ -f "${MERGED_PASS_VCF}" ]; then
    echo "  Final merged and PASS-filtered VCF: ${ABS_INPUT_VCF_DIR}/${MERGED_PASS_VCF}"
else
    echo "  PASS-filtered VCF was not successfully created."
fi

echo ""
echo "Summary Counts:"
if [ -f "${MERGED_RAW_VCF}" ]; then
    RAW_COUNT=$(bcftools view -H "${MERGED_RAW_VCF}" | grep -vc '^#' || echo "N/A")
    echo "  Number of variants in original merged file (${MERGED_RAW_VCF}): ${RAW_COUNT}"
fi
if [ -f "${MERGED_PASS_VCF}" ]; then
    PASS_COUNT=$(bcftools view -H "${MERGED_PASS_VCF}" | grep -vc '^#' || echo "N/A")
    echo "  Number of variants in PASS-filtered file (${MERGED_PASS_VCF}): ${PASS_COUNT}"
fi
echo "======================================================"
