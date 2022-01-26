
from importlib.metadata import files
from lib2to3.pgen2.token import STAR
from time import time

#samples currently hardcoded
SAMPLES = ["Chrtob3_D6_0mM_1_FD", "Chrtob3_D6_16mM_1_FD", "Chrtob3_D6_32mM_1_FD"]

#threads currently hardcoded
THREADS = 10

#expected outputs
rule all:
    input:
        expand("OD_summary_deinterleave/{sample}_deinterleave_summary.txt", sample = SAMPLES),
        expand("OD_summary_fastp/{sample}_fastp_summary.txt", sample = SAMPLES),
        expand("OD_sortmerna/{sample}_ready_for_star_fwd_0.fq.gz", sample = SAMPLES),
        "OD_star_input/combined_genome_assemblies.fasta",
        "OD_star_database/C_tobin_and_cohort_gff_annotated_GenBank.star_2.7.9a",
        expand("OD_star_alignment/star_alignment_{sample}_log", sample = SAMPLES),
        expand("OD_star_alignment/{sample}_aligned.sort.sam", sample = SAMPLES),
        expand("OD_htseq_counts/{sample}", sample = SAMPLES)
        

rule deinterleave:
    "Run BBmap reformat.sh"
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "OD_deinterleave/{sample}_read1.fq",
        "OD_deinterleave/{sample}_read2.fq",
        "OD_summary_deinterleave/{sample}_deinterleave_summary.txt"
    shell:
        """
	module load bbmap
        reformat.sh \
        in={input} \
        out1={output[0]} \
        out2={output[1]} \
        2>> {output[2]}
        """

rule fastp:
    "Run fastp"
    input:
        "OD_deinterleave/{sample}_read1.fq",
        "OD_deinterleave/{sample}_read2.fq",
        "OD_summary_deinterleave/{sample}_deinterleave_summary.txt"
    output:
        "OD_fastp/{sample}.report.json",
        "OD_fastp/{sample}.report.html",
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz",
        "OD_summary_fastp/{sample}_fastp_summary.txt"
    threads: THREADS
    shell:
        """
        fastp \
        -i {input[0]} \
        -I {input[1]} \
        -c \
        -w {threads} \
        --cut_front \
        --cut_front_window_size=1 \
        --cut_front_mean_quality=3 \
        --cut_tail \
        --cut_tail_window_size=1 \
        --cut_tail_mean_quality=3 \
        --cut_right \
        --cut_right_window_size=15 \
        --cut_right_mean_quality=15 \
        -l 35 \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -j {output[0]} \
        -h {output[1]} \
        -o {output[2]} \
        -O {output[3]} \
        1>> {output[4]}
        """

#download rRNA fasta files first.
rule sortmerna:
    "Remove rRNA"
    input:
        "OD_fastp/{sample}.out.R1.fq.gz",
        "OD_fastp/{sample}.out.R2.fq.gz",
        "data/sortmerna/rfam-5.8s-database-id98.fasta",
        "data/sortmerna/rfam-5s-database-id98.fasta",
        "data/sortmerna/silva-arc-16s-id95.fasta",
        "data/sortmerna/silva-arc-23s-id98.fasta",
        "data/sortmerna/silva-bac-16s-id90.fasta",
        "data/sortmerna/silva-bac-23s-id98.fasta",
        "data/sortmerna/silva-euk-18s-id95.fasta",
        "data/sortmerna/silva-euk-28s-id98.fasta"
    output:
        "OD_sortmerna/{sample}_ready_for_star_fwd_0.fq.gz",
        "OD_sortmerna/{sample}_ready_for_star_rev_0.fq.gz"
    threads: THREADS
    shell:
        """
        sortmerna \
        --ref {input[2]} \
        --ref {input[3]} \
        --ref {input[4]} \
        --ref {input[5]} \
        --ref {input[6]} \
        --ref {input[7]} \
        --ref {input[8]} \
        --ref {input[9]} \
        --reads {input[0]} \
        --reads {input[1]} \
        --workdir OD_sortmerna \
        --kvdb OD_sortmerna/{wildcards.sample}_ \
        --idx-dir OD_sortmerna/{wildcards.sample}_ \
        --readb OD_sortmerna/{wildcards.sample}_ \
        --threads {threads} \
        --aligned OD_sortmerna/{wildcards.sample}_rRNA_reads \
        --other OD_sortmerna/{wildcards.sample}_ready_for_star \
        --fastx \
        --paired_in \
        -out2
        """


rule star_input:
    "Create concatenated files for star database creation"
    input:
        "data/genome_assemblies",
        "data/genome_annotations"
    output:
        "OD_star_input/combined_genome_assemblies.fasta",
        "OD_star_input/combined_annotations.gtf"
    threads: THREADS
    shell:
        """
        cat {input[0]}/* | grep -v "==>" > {output[0]}
        cat {input[1]}/* > {output[1]}
        """

rule star_database:
    "Create STAR database"
    input:
        "OD_star_input/combined_genome_assemblies.fasta",
        "OD_star_input/combined_annotations.gtf"
    output:
        "OD_star_database/star_database_log",
        directory("OD_star_database/C_tobin_and_cohort_gff_annotated_GenBank.star_2.7.9a")
    threads: THREADS
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeSAindexNbases 12 \
        --genomeDir OD_star_database/C_tobin_and_cohort_gff_annotated_GenBank.star_2.7.9a \
        --genomeFastaFiles {input[0]} \
        --sjdbGTFfile {input[1]} 2> {output[0]}
        """


rule star_alignment:
    "Run STAR Alignment"
    input:
        "OD_sortmerna/{sample}_ready_for_star_fwd_0.fq.gz",
        "OD_sortmerna/{sample}_ready_for_star_rev_0.fq.gz",
        "OD_star_database/C_tobin_and_cohort_gff_annotated_GenBank.star_2.7.9a"
    output:
        "OD_star_alignment/star_alignment_{sample}_log",
        "OD_star_alignment/{sample}_Aligned.out.sam"
    threads: 20
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads \
        --outFilterMultimapNmax 3 \
        --outSAMunmapped Within KeepPairs \
        --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
        --readFilesCommand zcat \
        --readFilesIn {input[0]} {input[1]}  \
        --genomeDir {input[2]} \
        --outFileNamePrefix OD_star_alignment/{wildcards.sample}_ 2> {output[0]}
        """

rule samtools:
    "Run Sort with samtools"
    input:
       "OD_star_alignment/{sample}_Aligned.out.sam"
    output:
        "OD_star_alignment/{sample}_aligned.sort.sam"
    shell:
        """
        samtools sort -n -O sam -T {wildcards.sample}.sort -o {output} {input}
        """

# Htseq time

rule htseq:
    "Run HTseq_count as stranded"
    input:
       "OD_star_alignment/{sample}_aligned.sort.sam",
       "OD_star_input/combined_annotations.gtf"
    output:
        "OD_htseq_counts/{sample}",
        "OD_htseq_counts/OD_seq_counts_{sample}_log"
    threads: THREADS
    shell:
        """
        htseq-count --idattr=gene_id --stranded=reverse \
        {input[0]} \
        {input[1]} \
        1> {output[0]} 2> {output[1]}
        """


