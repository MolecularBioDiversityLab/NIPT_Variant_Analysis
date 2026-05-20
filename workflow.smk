# =============================================================================
# GATK Best Practices - Germline Short Variant Discovery Pipeline
# Per-sample QC → Trimming → Alignment → Variant Calling (GVCF)
# Joint genotyping steps are handled by supporting_scripts/
# =============================================================================

import os

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------
_raw_dir = os.path.join(config["working_dir"], config["pre_qc"]["input_dir"])
SAMPLES = [
    f.replace(".fq.gz", "")
    for f in os.listdir(_raw_dir)
    if f.endswith(".fq.gz")
]

if not SAMPLES:
    raise ValueError(f"No .fq.gz files found in {_raw_dir}. Check working_dir and pre_qc.input_dir in config.")

# ---------------------------------------------------------------------------
# Helper lambdas
# ---------------------------------------------------------------------------
def wd(*parts):
    """Join working_dir with path parts from config."""
    return os.path.join(config["working_dir"], *parts)

# ---------------------------------------------------------------------------
# rule all — final targets
# ---------------------------------------------------------------------------
rule all:
    input:
        # Pre-trim QC
        expand(wd(config["pre_qc"]["output_dir"], "{sample}_fastqc.html"),  sample=SAMPLES),
        expand(wd(config["pre_stats"]["output_dir"], "{sample}.stats"),     sample=SAMPLES),
        # Post-trim QC
        expand(wd(config["post_qc"]["output_dir"], "{sample}_fastqc.html"), sample=SAMPLES),
        expand(wd(config["post_stats"]["output_dir"], "{sample}.stats"),    sample=SAMPLES),
        # Alignment QC
        expand(wd(config["pre_bamstats"]["output_dir"],  "{sample}.flagstat"), sample=SAMPLES),
        expand(wd(config["post_bamstats"]["output_dir"], "{sample}.flagstat"), sample=SAMPLES),
        # Final per-sample GVCFs
        expand(wd(config["haplotypecaller"]["output_dir"], "{sample}.g.vcf.gz"), sample=SAMPLES),


# =============================================================================
# QUALITY CONTROL — raw reads
# =============================================================================

rule pre_qc:
    """FastQC on raw reads."""
    input:
        read = wd(config["pre_qc"]["input_dir"], "{sample}.fq.gz")
    output:
        html = wd(config["pre_qc"]["output_dir"], "{sample}_fastqc.html"),
        zip  = wd(config["pre_qc"]["output_dir"], "{sample}_fastqc.zip")
    params:
        outdir = wd(config["pre_qc"]["output_dir"])
    threads: 1
    log:
        wd("logs", "pre_qc", "{sample}.log")
    shell:
        """
        mkdir -p {params.outdir}
        {config[tools][fastqc]} {input.read} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}
        """

rule pre_stats:
    """seqkit stats on raw reads."""
    input:
        read = wd(config["pre_stats"]["input_dir"], "{sample}.fq.gz")
    output:
        stats = wd(config["pre_stats"]["output_dir"], "{sample}.stats")
    log:
        wd("logs", "pre_stats", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        {config[tools][seqkit]} stats -a {input.read} > {output.stats} 2> {log}
        """


# =============================================================================
# TRIMMING
# =============================================================================

rule trim:
    """fastp adapter/quality trimming."""
    input:
        read = wd(config["trim"]["input_dir"], "{sample}.fq.gz")
    output:
        trimmed = wd(config["trim"]["output_dir"], "{sample}.fq.gz"),
        html    = wd(config["trim"]["output_dir"], "{sample}_fastp.html"),
        json    = wd(config["trim"]["output_dir"], "{sample}_fastp.json")
    params:
        q = config["trim"]["qc_threshold"]
    threads: config["num_threads"]
    log:
        wd("logs", "trim", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.trimmed})
        {config[tools][fastp]} \
            -i {input.read} \
            -o {output.trimmed} \
            -q {params.q} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            2> {log}
        """


# =============================================================================
# QUALITY CONTROL — trimmed reads
# =============================================================================

rule post_qc:
    """FastQC on trimmed reads."""
    input:
        read = wd(config["post_qc"]["input_dir"], "{sample}.fq.gz")
    output:
        html = wd(config["post_qc"]["output_dir"], "{sample}_fastqc.html"),
        zip  = wd(config["post_qc"]["output_dir"], "{sample}_fastqc.zip")
    params:
        outdir = wd(config["post_qc"]["output_dir"])
    threads: 1
    log:
        wd("logs", "post_qc", "{sample}.log")
    shell:
        """
        mkdir -p {params.outdir}
        {config[tools][fastqc]} {input.read} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}
        """

rule post_stats:
    """seqkit stats on trimmed reads."""
    input:
        read = wd(config["post_stats"]["input_dir"], "{sample}.fq.gz")
    output:
        stats = wd(config["post_stats"]["output_dir"], "{sample}.stats")
    log:
        wd("logs", "post_stats", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        {config[tools][seqkit]} stats -a {input.read} > {output.stats} 2> {log}
        """


# =============================================================================
# ALIGNMENT
# =============================================================================

rule alignment:
    """BWA-MEM alignment → coordinate-sorted BAM."""
    input:
        read = wd(config["alignment"]["input_dir"], "{sample}.fq.gz")
    output:
        bam = temp(wd(config["alignment"]["output_dir"], "{sample}_sorted.bam"))
    params:
        ref = config["reference_genome"]
    threads: config["num_threads"]
    log:
        wd("logs", "alignment", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        {config[tools][bwa]} mem \
            -t {threads} \
            {params.ref} \
            {input.read} \
            2> {log} \
        | {config[tools][samtools]} view -bSh - \
        | {config[tools][samtools]} sort -@ {threads} -o {output.bam}
        """

rule add_rg:
    """Add read group tags required by GATK."""
    input:
        bam = wd(config["add_rg"]["input_dir"], "{sample}_sorted.bam")
    output:
        bam = temp(wd(config["add_rg"]["output_dir"], "{sample}_sorted_rg.bam"))
    params:
        prefix = config["add_rg"]["rg_prefix"]
    log:
        wd("logs", "add_rg", "{sample}.log")
    shell:
        """
        {config[tools][gatk]} AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -RGID {wildcards.sample} \
            -RGLB {wildcards.sample} \
            -RGPL ILLUMINA \
            -RGPU {params.prefix} \
            -RGSM {wildcards.sample} \
            --VALIDATION_STRINGENCY LENIENT \
            2> {log}
        """

rule index_bam:
    """Index the read-group BAM."""
    input:
        bam = wd(config["add_rg"]["output_dir"], "{sample}_sorted_rg.bam")
    output:
        bai = wd(config["add_rg"]["output_dir"], "{sample}_sorted_rg.bam.bai")
    threads: config["num_threads"]
    log:
        wd("logs", "index_bam", "{sample}.log")
    shell:
        """
        {config[tools][samtools]} index -@ {threads} {input.bam} 2> {log}
        """

rule pre_bamstats:
    """flagstat before duplicate marking."""
    input:
        bam = wd(config["pre_bamstats"]["input_dir"], "{sample}_sorted_rg.bam"),
        bai = wd(config["add_rg"]["output_dir"],      "{sample}_sorted_rg.bam.bai")
    output:
        stats = wd(config["pre_bamstats"]["output_dir"], "{sample}.flagstat")
    log:
        wd("logs", "pre_bamstats", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        {config[tools][samtools]} flagstat {input.bam} > {output.stats} 2> {log}
        """

rule markduplicates:
    """GATK MarkDuplicates — marks but does NOT remove duplicates (GATK best practice)."""
    input:
        bam = wd(config["markduplicates"]["input_dir"], "{sample}_sorted_rg.bam"),
        bai = wd(config["add_rg"]["output_dir"],        "{sample}_sorted_rg.bam.bai")
    output:
        bam     = wd(config["markduplicates"]["output_dir"], "{sample}_sorted_markdup.bam"),
        metrics = wd(config["markduplicates"]["output_dir"], "{sample}.metrics")
    log:
        wd("logs", "markduplicates", "{sample}.log")
    shell:
        """
        {config[tools][gatk]} MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            2> {log}
        """

rule post_bamstats:
    """flagstat after duplicate marking."""
    input:
        bam = wd(config["post_bamstats"]["input_dir"], "{sample}_sorted_markdup.bam")
    output:
        stats = wd(config["post_bamstats"]["output_dir"], "{sample}.flagstat")
    log:
        wd("logs", "post_bamstats", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        {config[tools][samtools]} flagstat {input.bam} > {output.stats} 2> {log}
        """


# =============================================================================
# VARIANT CALLING — per sample GVCF
# =============================================================================

rule haplotypecaller:
    """GATK HaplotypeCaller in GVCF mode."""
    input:
        bam = wd(config["haplotypecaller"]["input_dir"], "{sample}_sorted_markdup.bam"),
        ref = config["reference_genome"]
    output:
        gvcf = wd(config["haplotypecaller"]["output_dir"], "{sample}.g.vcf.gz")
    threads: config["num_threads"]
    log:
        wd("logs", "haplotypecaller", "{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.gvcf})
        {config[tools][gatk]} HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            2> {log}
        """
