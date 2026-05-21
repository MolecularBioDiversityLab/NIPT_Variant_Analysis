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
        # Joint genotyping + VQSR — final filtered VCF
        wd(config["vqsr"]["output_vcf"]),
        # Final Annotated VCF (SnpEff + ANNOVAR)
        wd(config["annotation"]["output_dir"], f"final_annotated.{config['annotation']['annovar']['buildver']}_multianno.vcf")


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


# =============================================================================
# JOINT GENOTYPING & VARIANT RECALIBRATION
# =============================================================================

rule make_sample_map:
    """Build a tab-separated sample_map from all per-sample GVCFs."""
    input:
        gvcfs = expand(
            wd(config["haplotypecaller"]["output_dir"], "{sample}.g.vcf.gz"),
            sample=SAMPLES
        )
    output:
        sample_map = wd(config["joint_genotyping"]["genomicsdb_workspace"], "sample_map.txt")
    run:
        import os
        os.makedirs(os.path.dirname(output.sample_map), exist_ok=True)
        with open(output.sample_map, "w") as fh:
            for gvcf in sorted(input.gvcfs):
                sample_name = os.path.basename(gvcf).replace(".g.vcf.gz", "")
                fh.write(f"{sample_name}\t{gvcf}\n")


rule genomicsdb_import:
    """Consolidate per-sample GVCFs into a GenomicsDB workspace."""
    input:
        sample_map = wd(config["joint_genotyping"]["genomicsdb_workspace"], "sample_map.txt")
    output:
        workspace = directory(wd(config["joint_genotyping"]["genomicsdb_workspace"], "db"))
    params:
        interval_list = config["joint_genotyping"]["interval_list"],
        batch_size    = config["joint_genotyping"]["batch_size"]
    log:
        wd("logs", "genomicsdb_import", "genomicsdb_import.log")
    shell:
        """
        {config[tools][gatk]} GenomicsDBImport \
            --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            --sample-name-map {input.sample_map} \
            --genomicsdb-workspace-path {output.workspace} \
            --batch-size {params.batch_size} \
            --genomicsdb-shared-posixfs-optimizations true \
            -L {params.interval_list} \
            2> {log}
        """


rule genotype_gvcfs:
    """Joint genotyping across all samples from the GenomicsDB workspace."""
    input:
        workspace = wd(config["joint_genotyping"]["genomicsdb_workspace"], "db"),
        ref       = config["reference_genome"]
    output:
        vcf = wd(config["joint_genotyping"]["output_vcf"])
    params:
        interval_list = config["joint_genotyping"]["interval_list"]
    log:
        wd("logs", "genotype_gvcfs", "genotype_gvcfs.log")
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {config[tools][gatk]} GenotypeGVCFs \
            --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            -R {input.ref} \
            -V gendb://{input.workspace} \
            -O {output.vcf} \
            -L {params.interval_list} \
            --genomicsdb-shared-posixfs-optimizations true \
            2> {log}
        """


rule variant_recalibration:
    """Build SNP and INDEL VQSR recalibration models (GATK best practice)."""
    input:
        vcf = wd(config["joint_genotyping"]["output_vcf"]),
        ref = config["reference_genome"]
    output:
        snp_recal    = wd(config["vqsr"]["output_prefix"] + ".snp.recal"),
        snp_tranches = wd(config["vqsr"]["output_prefix"] + ".snp.tranches"),
        snp_rscript  = wd(config["vqsr"]["output_prefix"] + ".snp.plots.R"),
        indel_recal    = wd(config["vqsr"]["output_prefix"] + ".indel.recal"),
        indel_tranches = wd(config["vqsr"]["output_prefix"] + ".indel.tranches"),
        indel_rscript  = wd(config["vqsr"]["output_prefix"] + ".indel.plots.R")
    params:
        hapmap = config["vqsr"]["snp"]["resources"]["hapmap"],
        omni   = config["vqsr"]["snp"]["resources"]["omni"],
        g1000  = config["vqsr"]["snp"]["resources"]["g1000"],
        dbsnp  = config["vqsr"]["snp"]["resources"]["dbsnp"],
        mills  = config["vqsr"]["indel"]["resources"]["mills"]
    log:
        snp   = wd("logs", "vqsr", "recalibration_snp.log"),
        indel = wd("logs", "vqsr", "recalibration_indel.log")
    shell:
        """
        mkdir -p $(dirname {output.snp_recal})
        mkdir -p $(dirname {log.snp})

        # SNP recalibration model
        {config[tools][gatk]} VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0  {params.hapmap} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0   {params.omni}   \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0  {params.g1000}  \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0   {params.dbsnp}  \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -O {output.snp_recal} \
            --tranches-file {output.snp_tranches} \
            --rscript-file  {output.snp_rscript} \
            2> {log.snp}

        # INDEL recalibration model
        {config[tools][gatk]} VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0  {params.mills}  \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0  {params.dbsnp}  \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL \
            -O {output.indel_recal} \
            --tranches-file {output.indel_tranches} \
            --rscript-file  {output.indel_rscript} \
            2> {log.indel}
        """


rule apply_vqsr:
    """Apply SNP then INDEL recalibration filters; emit the final VCF."""
    input:
        vcf          = wd(config["joint_genotyping"]["output_vcf"]),
        ref          = config["reference_genome"],
        snp_recal    = wd(config["vqsr"]["output_prefix"] + ".snp.recal"),
        snp_tranches = wd(config["vqsr"]["output_prefix"] + ".snp.tranches"),
        indel_recal    = wd(config["vqsr"]["output_prefix"] + ".indel.recal"),
        indel_tranches = wd(config["vqsr"]["output_prefix"] + ".indel.tranches")
    output:
        final_vcf = wd(config["vqsr"]["output_vcf"])
    params:
        snp_ts   = config["vqsr"]["snp"]["ts_filter_level"],
        indel_ts = config["vqsr"]["indel"]["ts_filter_level"],
        tmp_vcf  = lambda wc, output: output.final_vcf.replace(".vcf.gz", ".snp_recal.vcf.gz")
    log:
        snp   = wd("logs", "apply_vqsr", "apply_snp.log"),
        indel = wd("logs", "apply_vqsr", "apply_indel.log")
    shell:
        """
        mkdir -p $(dirname {output.final_vcf})
        mkdir -p $(dirname {log.snp})

        # Apply SNP VQSR
        {config[tools][gatk]} ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            -O {params.tmp_vcf} \
            --recal-file    {input.snp_recal} \
            --tranches-file {input.snp_tranches} \
            --ts-filter-level {params.snp_ts} \
            -mode SNP \
            2> {log.snp}

        # Apply INDEL VQSR
        {config[tools][gatk]} ApplyVQSR \
            -R {input.ref} \
            -V {params.tmp_vcf} \
            -O {output.final_vcf} \
            --recal-file    {input.indel_recal} \
            --tranches-file {input.indel_tranches} \
            --ts-filter-level {params.indel_ts} \
            -mode INDEL \
            2> {log.indel}

        # Remove intermediate SNP-only recalibrated VCF
        rm -f {params.tmp_vcf} {params.tmp_vcf}.tbi
        """


# =============================================================================
# VARIANT ANNOTATION (SnpEff & ANNOVAR)
# =============================================================================

rule snpeff_annotation:
    """Annotate variants functionally using SnpEff."""
    input:
        vcf = wd(config["vqsr"]["output_vcf"])
    output:
        vcf = wd(config["annotation"]["output_dir"], "final_filtered.snpeff.vcf")
    params:
        genome   = config["annotation"]["snpeff"]["genome_version"],
        data_dir = config["annotation"]["snpeff"]["data_dir"]
    log:
        wd("logs", "annotation", "snpeff.log")
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {config[tools][snpeff]} \
            -Xmx8g \
            -dataDir {params.data_dir} \
            {params.genome} \
            {input.vcf} > {output.vcf} 2> {log}
        """

rule annovar_annotation:
    """Further annotate the SnpEff output using ANNOVAR with dbSNP, ClinVar, and gnomAD."""
    input:
        vcf = wd(config["annotation"]["output_dir"], "final_filtered.snpeff.vcf")
    output:
        vcf = wd(config["annotation"]["output_dir"], f"final_annotated.{config['annotation']['annovar']['buildver']}_multianno.vcf"),
        txt = wd(config["annotation"]["output_dir"], f"final_annotated.{config['annotation']['annovar']['buildver']}_multianno.txt")
    params:
        humandb   = config["annotation"]["annovar"]["humandb_dir"],
        buildver  = config["annotation"]["annovar"]["buildver"],
        protocols = config["annotation"]["annovar"]["protocols"],
        ops       = config["annotation"]["annovar"]["operations"],
        prefix    = wd(config["annotation"]["output_dir"], "final_annotated")
    log:
        wd("logs", "annotation", "annovar.log")
    shell:
        """
        {config[tools][table_annovar]} {input.vcf} {params.humandb} \
            -buildver {params.buildver} \
            -out {params.prefix} \
            -remove \
            -protocol {params.protocols} \
            -operation {params.ops} \
            -nastring . \
            -vcfinput \
            2> {log}
        """
