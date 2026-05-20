#!/usr/bin/env bash
# =============================================================================
# supporting_scripts/02_genotype_gvcfs.sh
#
# Joint genotyping across all samples from a GenomicsDB workspace using
# GATK GenotypeGVCFs.
#
# Usage:
#   bash 02_genotype_gvcfs.sh \
#       <reference.fasta>   \
#       <genomicsdb_workspace> \
#       <output.vcf.gz>        \
#       <interval_list>
#
# Example:
#   bash 02_genotype_gvcfs.sh \
#       /ref/genome.fasta \
#       my_database \
#       genotyped.vcf.gz \
#       intervals.list
# =============================================================================
set -euo pipefail

REFERENCE="${1:?ERROR: reference FASTA required as first argument}"
WORKSPACE="${2:?ERROR: GenomicsDB workspace path required as second argument}"
OUTPUT_VCF="${3:?ERROR: output VCF path required as third argument}"
INTERVAL_LIST="${4:?ERROR: interval list required as fourth argument}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting GenotypeGVCFs"
echo "  Reference  : ${REFERENCE}"
echo "  Workspace  : ${WORKSPACE}"
echo "  Output     : ${OUTPUT_VCF}"
echo "  Intervals  : ${INTERVAL_LIST}"

mkdir -p "$(dirname "${OUTPUT_VCF}")"

gatk GenotypeGVCFs \
    --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R "${REFERENCE}" \
    -V "gendb://${WORKSPACE}" \
    -O "${OUTPUT_VCF}" \
    -L "${INTERVAL_LIST}" \
    --genomicsdb-shared-posixfs-optimizations true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] GenotypeGVCFs complete → ${OUTPUT_VCF}"
