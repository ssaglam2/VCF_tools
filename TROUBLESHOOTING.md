# Troubleshooting Guide for VCF Processing Pipeline

This guide provides solutions and debugging tips for common issues encountered while using the script.

## Table of Contents

1.  [Prerequisites Not Met](#1-prerequisites-not-met)
    *   [bcftools Not Found](#bcftools-not-found)
2.  [Script Execution Issues](#2-script-execution-issues)
    *   [Permission Denied](#permission-denied)
    *   [Bad Interpreter / Script Not Running](#bad-interpreter--script-not-running)
    *   [dos2unix Issues (Windows Line Endings)](#dos2unix-issues-windows-line-endings)
3.  [File Handling & Configuration Problems](#3-file-handling--configuration-problems)
    *   [No VCF Files Found or Processed](#no-vcf-files-found-or-processed)
    *   [Incorrect Sample Names in Merged VCF](#incorrect-sample-names-in-merged-vcf)
    *   [Intermediate Files Not Created as Expected](#intermediate-files-not-created-as-expected)
4.  [bcftools Command Errors](#4-bcftools-command-errors)
    *   [Errors During Indexing](#errors-during-indexing)
    *   [Errors During Reheadering](#errors-during-reheadering)
    *   [Errors During Merging](#errors-during-merging)
    *   [Errors During Filtering (Chromosome or PASS)](#errors-during-filtering-chromosome-or-pass)
5.  [Chromosome Filtering Specific Issues](#5-chromosome-filtering-specific-issues)
    *   [No Variants After Chromosome Filtering](#no-variants-after-chromosome-filtering)
6.  [General Debugging Tips](#6-general-debugging-tips)

---

## 1. Prerequisites Not Met

### bcftools Not Found

*   **Error Message:** `[FATAL ERROR] bcftools command not found...` or similar `bcftools: command not found`.
*   **Cause:** `bcftools` is either not installed or its installation directory is not in your system's `PATH` environment variable.
*   **Solution:**
    1.  **Verify Installation:** Type `bcftools --version` in your terminal. If it's not found, you need to install it.
    2.  **Install bcftools:** The recommended way is using Conda:
        ```bash
        conda create --name vcf_env -c bioconda bcftools
        conda activate vcf_env
        ```
        Alternatively, you can compile from source or use system package managers (e.g., `apt-get install bcftools` on Debian/Ubuntu), but Conda often manages dependencies better.
    3.  **Activate Conda Environment:** If installed via Conda, ensure the environment is activated *before* running the script: `conda activate [your_vcf_env_name]`.
    4.  **Check PATH:** If installed manually, ensure the directory containing the `bcftools` executable is added to your `PATH`.

## 2. Script Execution Issues

### Permission Denied

*   **Error Message:** `bash: ./run_vcf_processing.sh: Permission denied`
*   **Cause:** The script file does not have execute permissions.
*   **Solution:** Add execute permissions:
    ```bash
    chmod +x run_vcf_processing.sh
    ```

### Bad Interpreter / Script Not Running

*   **Error Message:** May vary, e.g., `bad interpreter: No such file or directory` or script tries to run with `sh` instead of `bash`, or commands within the script fail unexpectedly.
*   **Cause:**
    *   The shebang line (`#!/bin/bash`) at the top of the script might be incorrect or missing.
    *   The script might have Windows-style line endings (CRLF) instead of Unix-style (LF).
*   **Solution:**
    1.  **Check Shebang:** Ensure the very first line of the script is `#!/bin/bash`.
    2.  **Convert Line Endings:** Use `dos2unix` (recommended) or `sed` to convert line endings if you suspect this issue (especially if the script was edited on Windows):
        ```bash
        # Recommended
        dos2unix run_vcf_processing.sh

        # Alternative using sed
        # sed -i 's/\r$//' run_vcf_processing.sh
        ```

### dos2unix Issues (Windows Line Endings)
*   **Symptom:** Script fails with cryptic errors, commands behave strangely, or parts of lines are misinterpreted. `nano` might show `(Converted from DOS format)`.
*   **Cause:** Windows line endings (`\r\n`) instead of Unix (`\n`). The `\r` (carriage return) can break commands.
*   **Solution:** Run `dos2unix scriptname.sh` before executing. If `dos2unix` is not installed, use `sudo apt install dos2unix` (Debian/Ubuntu) or an equivalent for your system.

## 3. File Handling & Configuration Problems

### No VCF Files Found or Processed

*   **Symptom:** Script reports "Total candidates found (*.vcf.gz): 0" or "No VCF files were successfully prepared for merging."
*   **Cause & Solution:**
    1.  **Incorrect `INPUT_VCF_DIR`:** Verify that the `INPUT_VCF_DIR` variable in the script's configuration section points to the correct directory containing your `.vcf.gz` files. The script navigates to this directory.
    2.  **Files Not Ending in `.vcf.gz`:** The script specifically looks for files ending with `.vcf.gz`. Ensure your files have this exact extension (case-sensitive on Linux).
    3.  **`EXPECTED_FILENAME_SUBSTRING` Mismatch:**
        *   If `EXPECTED_FILENAME_SUBSTRING` is set in the configuration, the script will only process `.vcf.gz` files whose names *contain* this substring.
        *   Double-check that this substring accurately matches a part of all your target filenames. Pay attention to case sensitivity, hyphens, dots, etc.
        *   If you want to process *all* `.vcf.gz` files in the directory, set `EXPECTED_FILENAME_SUBSTRING=""` (empty string).
    4.  **No Files in Directory:** Ensure the `INPUT_VCF_DIR` actually contains the VCF files you intend to process. Use `ls -l` in that directory to check.

### Incorrect Sample Names in Merged VCF

*   **Symptom:** The sample names in the final merged VCF are not what you expected (e.g., they include too much of the filename or too little).
*   **Cause & Solution:** This is controlled by how the `NEW_SAMPLE_NAME_RAW` is derived using `basename` and the `ACTUAL_BASENAME_STRIP_SUFFIX` variable.
    1.  **Review `ACTUAL_BASENAME_STRIP_SUFFIX` Logic:**
        *   The script determines `ACTUAL_BASENAME_STRIP_SUFFIX` based on `BASENAME_STRIP_SUFFIX_OVERRIDE` and `EXPECTED_FILENAME_SUBSTRING`.
        *   If `BASENAME_STRIP_SUFFIX_OVERRIDE` is set, that value is used directly.
        *   Else, if `EXPECTED_FILENAME_SUBSTRING` is set, that is used as the suffix to strip.
        *   Else (if both are empty), `.vcf.gz` is stripped.
    2.  **Adjust Configuration:**
        *   The most common way to control this is by setting `EXPECTED_FILENAME_SUBSTRING` correctly. For example, if your files are `SAMPLEID.hard-filtered.vcf.gz` and you want `SAMPLEID`, set `EXPECTED_FILENAME_SUBSTRING=".hard-filtered.vcf.gz"`.
        *   If this doesn't work, use `BASENAME_STRIP_SUFFIX_OVERRIDE` to specify the exact suffix string you want `basename` to remove.
    3.  **Check Script Output:** The script prints the "Current internal sample" and the "New sanitized sample name" for each file. Review this output to see how names are being transformed.

### Intermediate Files Not Created as Expected

*   **Symptom:** The `REHEADERED_SUBDIR` is empty or missing expected `.renamed.vcf.gz` files.
*   **Cause & Solution:**
    1.  Usually a consequence of "No VCF Files Found or Processed" (see above). If no input files are processed, no intermediate files will be created.
    2.  Errors during the reheadering step for individual files (check script output for `ERROR: Failed to reheader...`).
    3.  Permissions issues preventing directory creation or file writing in `INPUT_VCF_DIR`.

## 4. bcftools Command Errors

The script logs which `bcftools` command is being run. If a `bcftools` command fails, it will typically print an error message.

### Errors During Indexing

*   **Symptom:** Messages like `ERROR: Failed to index ...`
*   **Cause:**
    *   Input VCF file is corrupted or not a valid VCF.
    *   Permissions issues preventing writing the `.tbi` index file.
    *   The VCF might not be bgzipped (though the script assumes `.vcf.gz`).
*   **Solution:**
    1.  Validate the problematic VCF file: `bcftools view your_file.vcf.gz > /dev/null`. If it errors, the VCF is likely problematic.
    2.  Check file permissions.

### Errors During Reheadering

*   **Symptom:** Messages like `ERROR: Failed to reheader ...`
*   **Cause:**
    *   Input VCF for reheadering is not indexed. (The script attempts to index it first).
    *   The temporary rename map file (`tmp_rename_map...txt`) is malformed or inaccessible.
    *   The new sample name being assigned is invalid (e.g., empty, or contains problematic characters not caught by sanitization, though unlikely).
*   **Solution:**
    1.  Check previous indexing steps for the input VCF.
    2.  Examine the script's logic for creating the rename map if this error is persistent.

### Errors During Merging

*   **Symptom:** `[ERROR] bcftools merge command failed.`
*   **Cause:**
    *   One or more input VCFs for merging are corrupted, not indexed, or have incompatible headers (e.g., different contig definitions if no reference FASTA is used by `bcftools merge`, though this script doesn't explicitly use one).
    *   Conflicting sample names (the script tries to avoid this by renaming).
    *   Incompatible VCF versions (rare).
*   **Solution:**
    1.  The script lists successfully prepared files. If only a few files are problematic, try removing them from the merge temporarily to isolate the issue.
    2.  Validate individual VCFs that were inputs to the merge.
    3.  Check `bcftools` error messages for more specific clues (e.g., "contig X not found in header").

### Errors During Filtering (Chromosome or PASS)

*   **Symptom:** `[ERROR] bcftools view for ... filtering failed.`
*   **Cause:**
    *   The input merged VCF for filtering is corrupted or not indexed. (The script attempts to index).
    *   The filter expression is invalid (for PASS filter, `'FILTER="PASS"'` is standard; for chromosome filter, region format is key).
*   **Solution:**
    1.  Ensure the VCF file being fed into the filter step exists and is valid.
    2.  For chromosome filtering, double-check the `CHROMS_TO_KEEP` format and names.

## 5. Chromosome Filtering Specific Issues

### No Variants After Chromosome Filtering

*   **Symptom:** The `*_autosomes_1_22.vcf.gz` (or similar) file is created but is empty or has far fewer variants than expected.
*   **Cause & Solution:**
    1.  **Chromosome Name Mismatch:** This is the most common cause. The names in `CHROMS_TO_KEEP` (e.g., `"1 2 3"`) **must exactly match** how chromosomes are named in your VCF header (e.g., contig lines might show `chr1`, `chr2`, `chr3`).
        *   **Action:** View your VCF header: `bcftools view -h your_merged_file.vcf.gz | grep '^##contig'`
        *   Adjust the `CHROMS_TO_KEEP` variable in the script to match your VCF's naming convention (e.g., change `"1"` to `"chr1"` if needed).
    2.  **No Variants on Those Chromosomes:** It's possible (though less likely for 1-22) that your merged VCF genuinely has no variants on the specified chromosomes that also pass other implicit filters (like being a variant call).

## 6. General Debugging Tips

*   **Run with `bash -x`:** Execute the script with `bash -x ./run_vcf_processing.sh`. This will print every command the script runs before it executes it, along with variable expansions. This is very verbose but extremely helpful for seeing exactly what the script is doing.
*   **Add `echo` Statements:** Insert `echo "DEBUG: Variable X = $X"` lines at various points in the script to check the values of variables or to see if certain code blocks are being reached.
*   **Test Commands Manually:** If a `bcftools` command in the script fails, copy that exact command (with filenames substituted) and run it directly in your terminal to see the full error message from `bcftools`.
*   **Check File Permissions:** Ensure you have read access to input files and write access to output directories.
*   **Check Disk Space:** Ensure you have sufficient disk space for intermediate and final VCF files, which can be large.
*   **Simplify:** If the full pipeline fails, try running only Stage 1, then only Stage 1 & 2, etc., to isolate where the problem occurs. Comment out later stages temporarily.
*   **Review Script Output Carefully:** The script tries to log what it's doing and any errors. Read these messages thoroughly.

---

If you encounter an issue not covered here, please try to gather as much information as possible (exact error messages, relevant configuration, `bcftools` version, sample of input filenames) when seeking help.
