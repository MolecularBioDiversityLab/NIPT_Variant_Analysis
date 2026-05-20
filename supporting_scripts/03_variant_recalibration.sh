#!/usr/bin/env bash
# =============================================================================
# supporting_scripts/03_variant_recalibration.sh
#
# GATK VariantRecalibrator — builds recalibration models for SNPs and INDELs
# separately (GATK best practice: run in two modes).
#
# Usage:
#   bash 03_variant_recalibration.sh \
#       <reference.fasta> \
#       <input.vcf.gz>    \
#       <output_prefix>   \
#       <hapmap.vcf.gz>   \
#       <omni.vcf.gz>     \
#       <1000G.vcf.gz>    \
#       <dbsnp.vcf.gz>    \
#       <mills.vcf.gz>
#
# Output files (per mode):
#   <output_prefix>.snp.recal
#   <output_prefix>.snp.tranches
#   <output_prefix>.snp.plots.R
#   <output_prefix>.indel.recal
#   <output_prefix>.indel.tranches
#   <output_prefix>.indel.plots.R
# =============================================================================
set -euo pipefail

REFERENCE="${1:?ERROR: reference FASTA required}"
INPUT_VCF="${2:?ERROR: input VCF required}"
OUT_PREFIX="${3:?ERROR: output prefix required}"
HAPMAP="${4:?ERROR: HapMap resource VCF required}"
OMNI="${5:?ERROR: 1000G Omni resource VCF required}"
G1000="${6:?ERROR: 1000G phase1 SNP resource VCF required}"
DBSNP="${7:?ERROR: dbSNP resource VCF required}"
MILLS="${8:?ERROR: Mills INDEL resource VCF required}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting VariantRecalibrator — SNP mode"

# ----- SNP recalibration ---------------------------------------------------
gatk VariantRecalibrator \
    -R "${REFERENCE}" \
    -V "${INPUT_VCF}" \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0   "${HAPMAP}" \
    --resource:omni,known=false,training=true,truth=false,prior=12.0    "${OMNI}"   \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0   "${G1000}"  \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0    "${DBSNP}"  \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O "${OUT_PREFIX}.snp.recal" \
    --tranches-file "${OUT_PREFIX}.snp.tranches" \
    --rscript-file  "${OUT_PREFIX}.snp.plots.R"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] VariantRecalibrator SNP done"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting VariantRecalibrator — INDEL mode"

# ----- INDEL recalibration -------------------------------------------------
gatk VariantRecalibrator \
    -R "${REFERENCE}" \
    -V "${INPUT_VCF}" \
    --resource:mills,known=false,training=true,truth=true,prior=12.0  "${MILLS}"  \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0  "${DBSNP}"  \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    -O "${OUT_PREFIX}.indel.recal" \
    --tranches-file "${OUT_PREFIX}.indel.tranches" \
    --rscript-file  "${OUT_PREFIX}.indel.plots.R"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] VariantRecalibrator INDEL done"
