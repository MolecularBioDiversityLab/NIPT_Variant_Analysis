# Short Variant Discovery Pipeline for NIPT Data

A Snakemake workflow implementing per-sample processing steps of the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) for Germline Short Variant Discovery. It handles raw sequencing data (FASTQ) and processes it through quality control, trimming, alignment, and variant calling to produce per-sample GVCF files.
---

## Overview

```
Raw FASTQ → QC → Trimming → Alignment → Duplicate Marking → HaplotypeCaller (GVCF)
                                                                        ↓
                                              GenomicsDBImport → GenotypeGVCFs → VQSR → Final VCF
```

The **Snakemake workflow** (`workflow.smk`) handles per-sample processing through GVCF generation.  
The **supporting scripts** handle joint genotyping and variant recalibration, which are typically run once after all samples have been processed.

---

## Repository Structure

```
.
├── workflow.smk                    # Main Snakemake workflow
├── config/
│   └── config.yaml                 # All user-configurable parameters
├── supporting_scripts/
│   ├── 01_genomicsdb_import.sh     # Consolidate GVCFs → GenomicsDB
│   ├── 02_genotype_gvcfs.sh        # Joint genotyping
│   ├── 03_variant_recalibration.sh # Build SNP + INDEL VQSR models
│   └── 04_apply_vqsr.sh            # Apply VQSR filters
└── README.md
```

---

## Requirements

| Tool | Version tested | Notes |
|------|---------------|-------|
| Snakemake | ≥ 7.0 | |
| GATK | ≥ 4.3 | |
| BWA | ≥ 0.7.17 | Reference must be BWA-indexed |
| samtools | ≥ 1.15 | |
| fastp | ≥ 0.23 | |
| FastQC | ≥ 0.11 | |
| seqkit | ≥ 2.0 | |

All tools are expected to be on `$PATH` (or full paths set in `config/config.yaml`).

---

## Quick Start

### 1. Prepare your reference genome

```bash
bwa index reference.fasta
samtools faidx reference.fasta
gatk CreateSequenceDictionary -R reference.fasta
```

### 2. Configure the pipeline

Edit `config/config.yaml` and set all fields marked `# <-- SET THIS`:

```yaml
working_dir: "/path/to/your/project"
reference_genome: "/path/to/reference.fasta"
```

Place raw FASTQ files (single-end, `.fq.gz`) in `<working_dir>/01.raw/00.fastq/`.

> **Note:** The pipeline currently assumes **single-end** reads. For paired-end support, the `trim` and `alignment` rules require modification.

### 3. Run the Snakemake workflow (per-sample)

```bash
# Dry-run first
snakemake -s workflow.smk --configfile config/config.yaml -n

# Run with 4 parallel jobs
snakemake -s workflow.smk --configfile config/config.yaml -j 4
```

### 4. Run joint genotyping (after all GVCFs are ready)

```bash
# Build a sample map (sample_name<TAB>path/to/sample.g.vcf.gz)
find <working_dir>/04.variant_calling/00.gvcf -name "*.g.vcf.gz" \
  | awk -F'/' '{print substr($NF,1,length($NF)-7)"\t"$0}' \
  > sample_map.txt

# Step 1: GenomicsDB import
bash supporting_scripts/01_genomicsdb_import.sh \
    sample_map.txt \
    <working_dir>/04.variant_calling/01.genomicsdb \
    /path/to/intervals.list

# Step 2: Joint genotyping
bash supporting_scripts/02_genotype_gvcfs.sh \
    /path/to/reference.fasta \
    <working_dir>/04.variant_calling/01.genomicsdb \
    <working_dir>/04.variant_calling/02.genotyped/genotyped.vcf.gz \
    /path/to/intervals.list

# Step 3: VQSR model building
bash supporting_scripts/03_variant_recalibration.sh \
    /path/to/reference.fasta \
    <working_dir>/04.variant_calling/02.genotyped/genotyped.vcf.gz \
    <working_dir>/04.variant_calling/03.vqsr/cohort \
    /path/to/hapmap.vcf.gz \
    /path/to/omni.vcf.gz \
    /path/to/1000G.vcf.gz \
    /path/to/dbsnp.vcf.gz \
    /path/to/mills.vcf.gz

# Step 4: Apply VQSR
bash supporting_scripts/04_apply_vqsr.sh \
    /path/to/reference.fasta \
    <working_dir>/04.variant_calling/02.genotyped/genotyped.vcf.gz \
    <working_dir>/04.variant_calling/03.vqsr/cohort \
    <working_dir>/04.variant_calling/04.final/final_filtered.vcf.gz
```

---

## Output Directory Layout

```
<working_dir>/
├── 01.raw/
│   ├── 00.fastq/        # Input raw reads
│   ├── 01.qc/           # FastQC reports (raw)
│   └── 02.stats/        # seqkit stats (raw)
├── 02.clean/
│   ├── 00.fastq/        # Trimmed reads + fastp reports
│   ├── 01.qc/           # FastQC reports (trimmed)
│   └── 02.stats/        # seqkit stats (trimmed)
├── 03.alignment/
│   ├── 00.bam/          # Sorted, RG-tagged, duplicate-marked BAMs + indices
│   ├── 01.pre_stats/    # flagstat before dedup
│   └── 02.post_stats/   # flagstat after dedup
├── 04.variant_calling/
│   ├── 00.gvcf/         # Per-sample GVCFs
│   ├── 01.genomicsdb/   # GenomicsDB workspace
│   ├── 02.genotyped/    # Joint-genotyped VCF
│   ├── 03.vqsr/         # Recalibration models
│   └── 04.final/        # Final filtered VCF
└── logs/                # Per-rule log files
```

---

## Notes

- **Duplicate marking:** Duplicates are **marked but not removed**, following GATK best practices. GATK tools are aware of the `DUPLICATE` flag and handle these reads appropriately.
- **Single-end vs paired-end:** This pipeline is configured for single-end reads. Paired-end adaptation requires changes to the `trim` and `alignment` rules.
- **VQSR cohort size:** VQSR requires a sufficiently large cohort (typically ≥ 30 whole-genome samples). For smaller cohorts, consider using `--filter-expression` hard filters instead.
