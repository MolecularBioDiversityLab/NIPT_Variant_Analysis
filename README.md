# GATK Short Variant Discovery Pipeline — NIPT Population Genetics

A Snakemake workflow implementing the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) for germline short variant discovery, adapted for **non-invasive prenatal testing (NIPT) cell-free DNA data**. The pipeline processes raw sequencing reads end-to-end — from FASTQ to a cohort-level, VQSR-filtered VCF — enabling downstream population genetics analyses such as allele frequency estimation, ancestry inference, and variant burden studies across NIPT cohorts.

---

## Overview

```
Raw FASTQ
   │
   ├── FastQC + seqkit stats (raw reads)
   │
   ├── fastp trimming
   │
   ├── FastQC + seqkit stats (trimmed reads)
   │
   ├── BWA-MEM alignment → coordinate-sorted BAM
   │
   ├── Add read groups (GATK requirement)
   │
   ├── flagstat (pre-dedup)
   │
   ├── GATK MarkDuplicates
   │
   ├── flagstat (post-dedup)
   │
   └── GATK HaplotypeCaller (per-sample GVCF)
            │
            └── [all samples complete]
                     │
                     ├── GenomicsDBImport   (consolidate GVCFs)
                     │
                     ├── GenotypeGVCFs      (joint genotyping)
                     │
                     ├── VariantRecalibrator (VQSR — SNP + INDEL models)
                     │
                     └── ApplyVQSR          → final_filtered.vcf.gz
```

The entire pipeline is managed by a single Snakemake workflow (`workflow.smk`). All per-sample and cohort-level steps are integrated — no manual hand-off between stages is required.

---

## NIPT-Specific Considerations

Cell-free DNA from maternal plasma presents characteristics that distinguish NIPT data from standard germline sequencing:

- **Low fetal fraction** — fetal DNA typically comprises 5–20% of total cfDNA. Variant callers may flag fetal-origin variants as low-allele-fraction somatic calls. Downstream filtering should account for expected heterozygous allele frequencies well below 0.5.
- **Fragment length bias** — fetal cfDNA fragments are shorter (~143 bp) than maternal cfDNA (~167 bp). Fragment-length-aware tools or filtering steps may be applied post-pipeline for fetal-specific analyses.
- **Single-end sequencing** — the pipeline is currently configured for single-end reads, which is common in high-throughput NIPT settings. Paired-end adaptation requires changes to the `trim` and `alignment` rules (see Notes).
- **Duplicate marking** — duplicates are **marked but not removed**, following GATK best practices. For cfDNA, duplicate rates can be high; metrics files from MarkDuplicates should be monitored per sample.
- **VQSR cohort size** — VQSR requires a sufficiently large cohort (≥ 30 whole-genome samples recommended). For smaller NIPT cohorts, consider hard filters using `--filter-expression` instead.

---

## Repository Structure

```
.
├── workflow.smk                    # Full Snakemake workflow (per-sample + joint genotyping)
├── config/
│   └── config.yaml                 # All user-configurable parameters
├── supporting_scripts/             # Standalone bash scripts (superseded by workflow.smk)
│   ├── 01_genomicsdb_import.sh
│   ├── 02_genotype_gvcfs.sh
│   ├── 03_variant_recalibration.sh
│   └── 04_apply_vqsr.sh
└── README.md
```

> The scripts in `supporting_scripts/` are retained for reference and manual re-runs, but all steps are now integrated into `workflow.smk`.

---

## Requirements

| Tool | Version tested | Notes |
|------|----------------|-------|
| Snakemake | ≥ 7.0 | |
| GATK | ≥ 4.3 | |
| BWA | ≥ 0.7.17 | Reference must be BWA-indexed |
| samtools | ≥ 1.15 | |
| fastp | ≥ 0.23 | |
| FastQC | ≥ 0.11 | |
| seqkit | ≥ 2.0 | |

All tools must be on `$PATH`, or full paths set in `config/config.yaml`.

---

## Quick Start

### 1. Prepare the reference genome

```bash
bwa index reference.fasta
samtools faidx reference.fasta
gatk CreateSequenceDictionary -R reference.fasta
```

### 2. Configure the pipeline

Edit `config/config.yaml`. At minimum, set all fields marked `# <-- SET THIS`:

```yaml
working_dir:      "/path/to/your/project"
reference_genome: "/path/to/reference.fasta"
```

Also set paths for the interval list and all VQSR resource VCFs (HapMap, 1000G Omni, 1000G phase1, dbSNP, Mills).

### 3. Stage input reads

Place raw FASTQ files (single-end, `.fq.gz`) in:

```
<working_dir>/01.raw/00.fastq/
```

Each file should be named `<sample_name>.fq.gz`. Sample names are discovered automatically.

### 4. Run the pipeline

```bash
# Dry-run first to verify the DAG
snakemake -s workflow.smk --configfile config/config.yaml -n

# Execute with 4 parallel jobs
snakemake -s workflow.smk --configfile config/config.yaml -j 4
```

A single run processes all samples through GVCF generation, then performs joint genotyping and VQSR to produce the final cohort VCF.

---

## Output Directory Layout

```
<working_dir>/
├── 01.raw/
│   ├── 00.fastq/              # Input raw reads
│   ├── 01.qc/                 # FastQC reports (raw)
│   └── 02.stats/              # seqkit stats (raw)
├── 02.clean/
│   ├── 00.fastq/              # Trimmed reads + fastp reports
│   ├── 01.qc/                 # FastQC reports (trimmed)
│   └── 02.stats/              # seqkit stats (trimmed)
├── 03.alignment/
│   ├── 00.bam/                # Sorted, RG-tagged, duplicate-marked BAMs + indices
│   ├── 01.pre_stats/          # flagstat before dedup
│   └── 02.post_stats/         # flagstat after dedup
├── 04.variant_calling/
│   ├── 00.gvcf/               # Per-sample GVCFs (.g.vcf.gz)
│   ├── 01.genomicsdb/         # GenomicsDB workspace + sample_map.txt
│   ├── 02.genotyped/          # Joint-genotyped cohort VCF
│   ├── 03.vqsr/               # Recalibration models (.recal, .tranches, .R plots)
│   └── 04.final/              # Final VQSR-filtered VCF ← primary output
└── logs/                      # Per-rule log files
```

---

## Workflow Rules

| Rule | Description |
|------|-------------|
| `pre_qc` | FastQC on raw reads |
| `pre_stats` | seqkit stats on raw reads |
| `trim` | fastp adapter/quality trimming |
| `post_qc` | FastQC on trimmed reads |
| `post_stats` | seqkit stats on trimmed reads |
| `alignment` | BWA-MEM → coordinate-sorted BAM |
| `add_rg` | GATK AddOrReplaceReadGroups |
| `index_bam` | samtools index |
| `pre_bamstats` | samtools flagstat (before MarkDuplicates) |
| `markduplicates` | GATK MarkDuplicates |
| `post_bamstats` | samtools flagstat (after MarkDuplicates) |
| `haplotypecaller` | GATK HaplotypeCaller in GVCF mode |
| `make_sample_map` | Auto-generate sample_map.txt from all GVCFs |
| `genomicsdb_import` | GATK GenomicsDBImport |
| `genotype_gvcfs` | GATK GenotypeGVCFs (joint genotyping) |
| `variant_recalibration` | GATK VariantRecalibrator (SNP + INDEL models) |
| `apply_vqsr` | GATK ApplyVQSR → final_filtered.vcf.gz |

---

## Downstream Population Genetics

The final VCF (`04.variant_calling/04.final/final_filtered.vcf.gz`) is intended as input for population-level analyses. Typical next steps include:

- **Allele frequency estimation** — cohort-level AF/AC/AN annotations are emitted by GenotypeGVCFs and can be used directly for frequency spectrum analyses.
- **Ancestry inference** — project cohort variants onto reference panels (e.g., 1000 Genomes, gnomAD) using tools such as PLINK, KING, or PC-AiR for principal component analysis and ancestry assignment.
- **Hardy–Weinberg and population stratification** — standard QC steps with PLINK or hail before association or burden testing.
- **Fetal fraction-aware filtering** — for NIPT-specific analyses, variants with allele balance consistent with fetal-fraction heterozygosity may be flagged or stratified separately.

---

## Notes

- **Single-end vs paired-end:** The pipeline is configured for single-end reads. Paired-end adaptation requires changes to the `trim` and `alignment` rules (add `-I input_R2` / second input fastq and enable paired mode in fastp/BWA).
- **Duplicate marking:** Duplicates are marked but not removed, per GATK best practices. GATK tools are aware of the `DUPLICATE` flag.
- **VQSR cohort size:** VQSR requires a large cohort (≥ 30 WGS samples). For small cohorts, use hard filters via `gatk VariantFiltration` instead.
- **GenomicsDB workspace:** The workspace path in `config.yaml` (`joint_genotyping.genomicsdb_workspace`) must not already exist when the pipeline runs for the first time. Delete or rename it before re-running `genomicsdb_import`.
- **Intermediate files:** The SNP-only intermediate VCF produced during `apply_vqsr` is automatically deleted after INDEL VQSR is applied.
