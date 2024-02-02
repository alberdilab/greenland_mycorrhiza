### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fq.gz", "")
            for fn in glob(f"resources/reads/*_1.fq.gz")]

print("Detected the following samples:")
print(SAMPLE)

################################################################################
### Setup the desired outputs
rule all:
    input:
        "results/coverm/counts.tsv"
        
################################################################################
### Filter reads with fastp
rule fastp:
    input:
        r1 = "resources/reads/{sample}_1.fq.gz",
        r2 = "resources/reads/{sample}_2.fq.gz"
    output:
        r1 = temp("results/fastp/{sample}_1.fq.gz"),
        r2 = temp("results/fastp/{sample}_2.fq.gz"),
        fastp_html = "results/fastp/{sample}.html",
        fastp_json = "results/fastp/{sample}.json"
    conda:
        "environment.yaml"
    params:
        jobname = "fastp_{sample}",
    threads:
        10
    resources:
        mem_gb=24,
        time='01:00:00'
    log:
        "logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --adapter_sequence_r2 CTGTCTCTTATACACATCT \
        &> {log}
        """

################################################################################
## Index Unite database:
rule index_unite:
    input:
        "resources/database/sh_general_release_dynamic_25.07.2023.fasta"
    output:
        multiext(
            "resources/database/sh_general_release_dynamic_25.07.2023",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2")
    conda:
        "environment.yaml"
    params:
        jobname = "mags_index",
        database = "resources/database/sh_general_release_dynamic_25.07.2023"
    threads:
        24
    resources:
        mem_gb=24,
        time='02:00:00'
    log:
        "logs/unite_index.log"
    message:
        "Indexing unite database with Bowtie2"
    shell:
        """
        # Index MAG gene catalogue
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {params.database} \
        &> {log}


        """
        
################################################################################
### Map non-host reads to DRAM genes files using Bowtie2
rule bowtie2_unite_mapping:
    input:
        r1 = "results/fastp/{sample}_1.fq.gz",
        r2 = "results/fastp/{sample}_2.fq.gz",
        index=multiext(
            "resources/database/sh_general_release_dynamic_25.07.2023",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2")
    output:
        bam = "results/bowtie/{sample}.bam"
    params:
        jobname = "mags_{sample}",
        database = "resources/database/sh_general_release_dynamic_25.07.2023"
    conda:
        "environment.yaml"
    threads:
        20
    resources:
        mem_gb=24,
        time='02:00:00'
    log:
        "logs/{sample}_bowtie.log"
    message:
        "Mapping {wildcards.sample} to unite using Bowtie2"
    shell:
        """
        # Map reads to MAGs using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.database} \
            -1 {input.r1} \
            -2 {input.r2} \
            --seed 1337 \
        | samtools sort -@ {threads} -o {output.bam} \
        &> {log}
        """

################################################################################
### Calculate the number of reads that mapped to MAG catalogue genes with CoverM
rule coverM_counts:
    input:
        expand("results/bowtie/{sample}.bam", sample=SAMPLE),
    output:
        gene_counts = "results/coverm/counts.tsv",
    params:
        jobname = "coverm"
    conda:
        "environment.yaml"
    threads:
        24
    resources:
        mem_gb=24,
        time='01:00:00'
    message:
        "Calculating MAG gene mapping rate using CoverM"
    shell:
        """
        coverm contig \
            -b {input} \
            -m count \
            -t {threads} \
            --proper-pairs-only \
            > {output.gene_counts}
        """
