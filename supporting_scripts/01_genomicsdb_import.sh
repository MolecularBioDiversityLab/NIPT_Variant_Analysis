#!/usr/bin/env bash
# =============================================================================
# supporting_scripts/01_genomicsdb_import.sh
#
# Consolidate per-sample GVCFs into a GenomicsDB workspace using
# GATK GenomicsDBImport.
#
# Usage:
#   bash 01_genomicsdb_import.sh \
#       <sample_map>     \   # tab-separated: sample_name<TAB>path/to/sample.g.vcf.gz
#       <workspace>      \   # output GenomicsDB directory (must NOT already exist)
#       <interval_list>  \   # .interval_list or .bed file
#       [batch_size]         # optional, default 50
#
# Example:
#   bash 01_genomicsdb_import.sh sample_map.txt my_database intervals.list 50
# =============================================================================
set -euo pipefail

SAMPLE_MAP="${1:?ERROR: sample_map file required as first argument}"
WORKSPACE="${2:?ERROR: GenomicsDB workspace path required as second argument}"
INTERVAL_LIST="${3:?ERROR: interval list required as third argument}"
BATCH_SIZE="${4:-50}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting GenomicsDBImport"
echo "  Sample map   : ${SAMPLE_MAP}"
echo "  Workspace    : ${WORKSPACE}"
echo "  Intervals    : ${INTERVAL_LIST}"
echo "  Batch size   : ${BATCH_SIZE}"

gatk GenomicsDBImport \
    --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${WORKSPACE}" \
    --batch-size "${BATCH_SIZE}" \
    --genomicsdb-shared-posixfs-optimizations true \
    -L "${INTERVAL_LIST}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] GenomicsDBImport complete → ${WORKSPACE}"
