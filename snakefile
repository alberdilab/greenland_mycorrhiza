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
        "results/coverm/gene_counts.tsv"
        
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
## Index host genomes:
rule star_index:
     input:
            genome="resources/reference/host/GCF_003160815.1_VulVul2.2_genomic.fna",
            annotation="resources/reference/host/GCF_003160815.1_VulVul2.2_genomic.gtf"
     output:
            folder=directory("resources/reference/host/index")
     conda:
         "environment.yaml"
     params:
        jobname = "star_index",
     threads:
         24
     resources:
         mem_gb=96,
         time='02:00:00'
     log:
         "logs/star_index.log"
     message:
         "Indexing host genome"
     shell:
         """
        mkdir {output}
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.folder} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation} \
        2> {log} 1>&2
         """

################################################################################
### Map to host reference genome using STAR
rule star_mapping:
    input:
        r1 = "results/fastp/{sample}_1.fq.gz",
        r2 = "results/fastp/{sample}_2.fq.gz",
        index = "resources/reference/host/index"
    output:
        r1 = "results/star/{sample}_1.fq.gz",
        r2 = "results/star/{sample}_2.fq.gz",
        host_bam = "results/star/{sample}_host.bam"
    params:
        jobname = "star_{sample}",
        r1 = "results/star/{sample}_1.fq.",
        r2 = "results/star/{sample}_2.fq",
        gene_counts = "results/star/{sample}_read_counts.tsv",
        sj = "results/star/{sample}_SJ.tsv"
    conda:
        "environment.yaml"
    threads:
        24
    resources:
        mem_gb=24,
        time='02:00:00'
    log:
        "logs/{sample}_star.log"
    message:
        "Mapping {wildcards.sample} to host genome using STAR"
    shell:
        """
        # Map reads to host genome using STAR
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {wildcards.sample} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
        &> {log}

        # Rename files
        mv {wildcards.sample}Aligned.out.bam {output.host_bam}
        mv {wildcards.sample}ReadsPerGene.out.tab {params.gene_counts}
        mv {wildcards.sample}SJ.out.tab {params.sj}
        mv {wildcards.sample}Unmapped.out.mate1 {params.r1}
        mv {wildcards.sample}Unmapped.out.mate2 {params.r2}

        # Compress non-host reads
        pigz \
            -p {threads} \
            {output.r1}

        pigz \
            -p {threads} \
            {output.r1}
        """

################################################################################
### Functionally annotate MAGs with DRAM
rule dram:
    input:
        "resources/reference/microbiome/MAGs.fa"
    output:
        annotations = "results/dram/MAGs_annotations.tsv.gz",
        genes = "results/dram/MAGs_genes.fna.gz",
        genesfaa = "results/dram/MAGs_genes.faa.gz",
        genesgff = "results/dram/MAGs_genes.gff.gz",
        scaffolds = "results/dram/MAGs_scaffolds.fna.gz",
        gbk = "results/dram/MAGs.gbk.gz",
        distillate = directory("results/dram/MAGs_distillate")
    params:
        jobname = "dram",
        outdir = "results/dram/MAGs_annotate",
        mainout = "results/dram",
        trnas = "results/dram/MAGs_trnas.tsv.gz",
        rrnas = "results/dram/MAGs_rrnas.tsv.gz",
    threads:
        2
    resources:
        mem_gb=128,
        time='04:00:00'
    log:
        "logs/dram.log"
    message:
        "Using DRAM to functionally annotate mags"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules

        rm {params.outdir}

        DRAM.py annotate \
            -i {input} \
            -o {params.outdir} \
            --threads {threads} \
            --min_contig_size 1500 

        if test ! -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        echo "neither trnas nor rrnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            -o {output.distillate}

        elif test -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        echo "both trnas and rrnas are available"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}
        mv {params.outdir}/trnas.tsv {params.trnas}

        elif test -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        echo "only trnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/trnas.tsv {params.trnas}

        elif test ! -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        echo "only rrnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}

        fi    

        pigz -p {threads} {params.outdir}/*.tsv
        pigz -p {threads} {params.outdir}/*.fna
        pigz -p {threads} {params.outdir}/*.faa
        pigz -p {threads} {params.outdir}/*.gff
        pigz -p {threads} {params.outdir}/genbank/*
        pigz -p {threads} {output.distillate}/*

        mv {params.outdir}/annotations.tsv.gz {output.annotations}
        mv {params.outdir}/scaffolds.fna.gz {output.scaffolds}
        mv {params.outdir}/genes.fna.gz {output.genes}
        mv {params.outdir}/*.faa.gz {output.genesfaa}
        mv {params.outdir}/*.gff.gz {output.genesgff}
        mv {params.outdir}/genbank/* {output.gbk}

        """

################################################################################
## Index MAGs:
rule index_mags:
    input:
        "results/dram/MAGs_genes.fna.gz"
    output:
        "results/dram/MAGs_genes"
    conda:
        "environment.yaml"
    params:
        jobname = "mags_index",
    threads:
        24
    resources:
        mem_gb=24,
        time='02:00:00'
    log:
        "logs/MAG_genes_index.log"
    message:
        "Indexing MAG catalogue with Bowtie2"
    shell:
        """
        # Index MAG gene catalogue
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {output} \
        &> {log}
        """
        
################################################################################
### Map non-host reads to DRAM genes files using Bowtie2
rule bowtie2_MAG_mapping:
    input:
        r1 = "results/star/{sample}_1.fq.gz",
        r2 = "results/star/{sample}_2.fq.gz",
        bt2_index = "results/dram/MAGs_genes"
    output:
        bam = "results/bowtie/{sample}.bam"
    params:
        jobname = "mags_{sample}",
        MAG_genes = "resources/reference/microbiome/MAG_genes.fna.gz"
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
        "Mapping {wildcards.sample} to MAG genes using Bowtie2"
    shell:
        """
        # Map reads to MAGs using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.MAG_genes} \
            -1 {input.r1} \
            -2 {input.r2} \
            --seed 1337 \
        | samtools sort -@ {threads} -o {output.bam} \
        &> {log}
        """

################################################################################
### Calculate the number of reads that mapped to MAG catalogue genes with CoverM
rule coverM_MAG_genes:
    input:
        expand("results/bowtie/{sample}.bam", sample=SAMPLE),
    output:
        gene_counts = "results/coverm/gene_counts.tsv",
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
