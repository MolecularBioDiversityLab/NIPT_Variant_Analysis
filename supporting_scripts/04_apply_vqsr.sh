#!/usr/bin/env bash
# =============================================================================
# supporting_scripts/04_apply_vqsr.sh
#
# GATK ApplyVQSR — applies SNP and INDEL recalibration models sequentially
# (GATK best practice: apply SNP first, then INDEL on the SNP-filtered output).
#
# Usage:
#   bash 04_apply_vqsr.sh \
#       <reference.fasta>  \
#       <input.vcf.gz>     \
#       <out_prefix>       \   # prefix for recal/tranches files from step 03
#       <final_output.vcf.gz>  \
#       [snp_ts_filter]    \   # optional, default 99.5
#       [indel_ts_filter]      # optional, default 99.0
#
# Example:
#   bash 04_apply_vqsr.sh \
#       /ref/genome.fasta \
#       genotyped.vcf.gz \
#       vqsr/cohort \
#       final_filtered.vcf.gz
# =============================================================================
set -euo pipefail

REFERENCE="${1:?ERROR: reference FASTA required}"
INPUT_VCF="${2:?ERROR: input VCF required}"
OUT_PREFIX="${3:?ERROR: recal/tranches prefix required (same as used in step 03)}"
FINAL_OUTPUT="${4:?ERROR: final output VCF path required}"
SNP_TS="${5:-99.5}"
INDEL_TS="${6:-99.0}"

INTERMEDIATE_VCF="${FINAL_OUTPUT%.vcf.gz}.snp_recal.vcf.gz"

mkdir -p "$(dirname "${FINAL_OUTPUT}")"

# ----- Apply SNP recalibration ---------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying SNP VQSR (ts_filter_level=${SNP_TS})"

gatk ApplyVQSR \
    -R "${REFERENCE}" \
    -V "${INPUT_VCF}" \
    -O "${INTERMEDIATE_VCF}" \
    --recal-file    "${OUT_PREFIX}.snp.recal" \
    --tranches-file "${OUT_PREFIX}.snp.tranches" \
    --ts-filter-level "${SNP_TS}" \
    -mode SNP

echo "[$(date '+%Y-%m-%d %H:%M:%S')] SNP VQSR applied → ${INTERMEDIATE_VCF}"

# ----- Apply INDEL recalibration -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying INDEL VQSR (ts_filter_level=${INDEL_TS})"

gatk ApplyVQSR \
    -R "${REFERENCE}" \
    -V "${INTERMEDIATE_VCF}" \
    -O "${FINAL_OUTPUT}" \
    --recal-file    "${OUT_PREFIX}.indel.recal" \
    --tranches-file "${OUT_PREFIX}.indel.tranches" \
    --ts-filter-level "${INDEL_TS}" \
    -mode INDEL

echo "[$(date '+%Y-%m-%d %H:%M:%S')] INDEL VQSR applied → ${FINAL_OUTPUT}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline complete. Final VCF: ${FINAL_OUTPUT}"

# Clean up intermediate file
rm -f "${INTERMEDIATE_VCF}" "${INTERMEDIATE_VCF}.tbi"
