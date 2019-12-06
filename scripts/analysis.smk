#!/usr/bin/env snakemake
import os
import glob

configfile: "config.yaml"

reference = config["REFERENCE"]
index = config["REFERENCE"]
datafolder = config["datafolder"]
projectdir = config["projectdir"]

# get files from directories in datafolder path
directories_in_datafolder = [directory for directory in glob.glob(datafolder + "/SRR*")]

rnaseq=[]
hic=[]
dnase=[]
chipseq=[]
bisulfite=[]

datafolder_structure = dict()
# print(directories_in_datafolder)

for directory in directories_in_datafolder:
    SRR=os.path.basename(directory).split("--")[0]
    datafolder_structure[SRR]={"directory":directory}
    for root, subdir, files in os.walk(directory):
        for filename in files:
            if "_1.fastq" in filename:
                datafolder_structure[SRR].update({"R1": filename})
            elif "_2.fastq" in filename:
                datafolder_structure[SRR].update({"R2": filename})
    if "RNA-Seq" in directory:
        rnaseq.append(SRR)
    elif "ChIP-Seq" in directory:
        chipseq.append(SRR)
    elif "DNase-Hypersensitivity" in directory:
        dnase.append(SRR)
    elif "HiC" in directory:
        hic.append(SRR)
    elif "Bisulfite-Seq" in directory:
        bisulfite.append(SRR)
    else:
        pass

# print(datafolder_structure)
# get the folders and files in a dictionary datatype
#
# directory, SRA = glob_wildcards(datafolder + "/{dir}/{srr}_1.fastq.gz")
#
# for x in zip(directory, SRA):
#     print(x)

rule minimap2_index:
    input: config['REFERENCE']
    output: config['REFERENCE'] + ".minimap2.index"
    message: "Minimap2 Indexing the reference"
    benchmark: "{projectdir}/results/benchmark/minimap2_reference_indexing.txt".format(projectdir=projectdir)
    threads: 4
    shell: "minimap2 -k 21 -d {output} {input}"

rule minimap2_align:
    input:
        R1=lambda wildcards: datafolder_structure[wildcards.srr]["directory"] + "/" + datafolder_structure[wildcards.srr]["R1"],
        referenceidx=config['REFERENCE'] + ".minimap2.index"

    output: "{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam"
    message: "Minimap2 aligning reads"

    run:
        if os.path.exists(input.R1.replace("_1.fastq", "_2.fastq")):
            R2=input.R1.replace("_1.fastq", "_2.fastq")

            shell("minimap2 -t 4 -a -R '@RG\tID={wildcards.srr}\tSM:{wildcards.srr}' -x sr  {input.referenceidx} {input.R1} {R2} | samtools view -b | samtools sort -o {output}" )
        else:
            shell("minimap2 -t 4 -a -R '@RG\tID={wildcards.srr}\tSM:{wildcards.srr}' -x sr  {input.referenceidx} {input.R1} | samtools view -b | samtools sort -o {output}" )


rule bowtie2_build_index:
    input: config['REFERENCE']
    output:
        config['REFERENCE']+".1.bt2",
        config['REFERENCE']+".2.bt2",
        config['REFERENCE']+".3.bt2",
        config['REFERENCE']+".4.bt2",
        config['REFERENCE']+".rev.1.bt2",
        config['REFERENCE']+".rev.2.bt2"

    message: "Building index for the reference sequence"
    threads: 4
    shell: "bowtie2-build -f {input} {input}"

# rule rnaseq_align_paired:
# 	"""
# 	trinity align reads to reference
# 	"""
# 	input:
# 		R1=lambda wildcards: datafolder_structure[wildcards.srr]["directory"] + "/" + datafolder_structure[wildcards.srr]["R1"]
# 	output: "{projectdir}/results/{srr}_rnaseq_aligned.bam"
# 	run:
#         # look for R2 file, if exists it is paired data
#         if os.path.exists(input.R1.replace("_1.fastq", "_2.fastq")):
#             R2=input.R1.replace("_1.fastq", "_2.fastq")
#             shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -1 {input.R1} -2 {R2} -S {output}" )
#         else:
#             shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -U {input.R1} -S {output}")
#########################################################################################################################
rule tophat_reads_align:
    """
    Align RNAseq data using tophat
    """
    input:
        R1=lambda wildcards: datafolder_structure[wildcards.srr]["directory"] + "/" + datafolder_structure[wildcards.srr]["R1"]
        # R2=lambda wildcards: datafolder_structure[wildcards.srr]["R2"]
    output: "{projectdir}/results/{datatype}/{srr}/accepted_hits.bam"
    threads: 4
    run:
        # look for R2 file, if exists it is paired data
        if os.path.exists(input.R1.replace("_1.fastq", "_2.fastq")):
            R2=input.R1.replace("_1.fastq", "_2.fastq")
            shell("tophat2 --rg-id {wildcards.srr} --rg-sample {wildcards.srr} --output-dir {projectdir}/results/{datatype}/{wildcards.srr} {reference} {input.R1}  {R2} " )
        else:
            shell("tophat2 --rg-id {wildcards.srr} --rg-sample {wildcards.srr} --output-dir {projectdir}/results/{datatype}/{wildcards.srr} {reference} {input.R1}")

rule tophat_reads_align_merge:
    """
    Merge bam files to a sorted bam
    """
    input: expand("{projectdir}/results/{datatype}/{srr}/accepted_hits.bam", projectdir=projectdir, srr=rnaseq, datatype="rnaseq")
    output: "{projectdir}/results/{datatype}/rnaseq_reads_bowtie_aligned_merged_accepted_hits.bam"
    threads: 4
    message: "merging rna seq aligned bam files"
    shell: "samtools merge -r -l 4 {output} {input} &&  samtools index {output}"


rule run_rnaseq_bowtie:
    """
    run rnaseq alignment and merge to a sorted bam file
    """
    # input: expand("{projectdir}/results/{srr}/accepted_hits.bam", projectdir=projectdir, srr=rnaseq)
    # input:expand("{projectdir}/results/{datatype}/{srr}/accepted_hits.bam", projectdir=projectdir, datatype="rnaseq", srr=rnaseq)
    input: expand("{projectdir}/results/{datatype}/rnaseq_reads_bowtie_aligned_merged_accepted_hits.bam", datatype="rnaseq", projectdir=projectdir)

rule run_rnaseq_minimap2:
    input: expand("{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam", projectdir=projectdir, datatype="rnaseq", srr=rnaseq)
    # output: expand("{projectdir}/results/{datatype}/rnaseq_reads_minimap_aligned_merged_accepted_hits.bam", datatype="rnaseq", projectdir=projectdir)
    # benchmark: "{projectdir}/results/benchmark/minimap2_align_mergedbam.txt".format(projectdir=projectdir)
    # message: "Merging bam files aligned with minimap2"
    # threads: 4
    # shell: "samtools merge -l 5 {output} {input} && samtools index {output}"

######################################################################################################################################


rule bowtie2_reads_align:
    input:
        R1=lambda wildcards: datafolder_structure[wildcards.srr]["directory"] + "/" + datafolder_structure[wildcards.srr]["R1"]
        # R2=lambda wildcards: datafolder_structure[wildcards.srr]["R2"]
    output:
        "{projectdir}/results/{datatype}/{srr}_bowtie_aligned.bam"
    threads: 4
    run:
        # look for R2 file, if exists it is paired data
        if os.path.exists(input.R1.replace("_1.fastq", "_2.fastq")):
            R2=input.R1.replace("_1.fastq", "_2.fastq")
            shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -1 {input.R1} -2 {R2} | samtools view -b | samtools sort -o {output} && samtools index {output}")
        else:
            shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -U {input.R1} | samtools view -b | samtools sort -o {output} && samtools index {output}")

#rule bowtie2_reads_align_merge:
#
# rule bismark_reads_align:
#     input:
#         R1=lambda wildcards: datafolder_structure[wildcards.srr]["directory"] + "/" + datafolder_structure[wildcards.srr]["R1"]
#         # R2=lambda wildcards: datafolder_structure[wildcards.srr]["R2"]
#     output:
#         "{projectdir}/results/{datatype}/bisulfite_{srr}_bowtie_aligned.bam"
#     threads: 4
#     run:
#         # look for R2 file, if exists it is paired data
#         if os.path.exists(input.R1.replace("_1.fastq", "_2.fastq")):
#             R2=input.R1.replace("_1.fastq", "_2.fastq")
#             shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -1 {input.R1} -2 {R2} | samtools view -b | samtools sort -o {output}" )
#         else:
#             shell("bowtie2 --no-unal --no-discordant --rg-id {wildcards.srr} -x {index} -U {input.R1} | samtools view -b | samtools sort -o {output}")
#
#
# rule bowtie2_unpaired_align:
#     input:"{dir}/{srr}.fastq.gz"
#     output: temp("{projectdir}/results/{srr}_bowtie_align.sam")
#     message: "Aligning single end data with bowtie2 "
#     threads: 2
#     shell: "bowtie2 --rg-id {wildcards.srr} -x {index} -U {input} -S {output}"



rule run_dnase_bowtie:
    input: expand("{projectdir}/results/{datatype}/{srr}_bowtie_aligned.bam", projectdir=projectdir, datatype="dnase", srr=dnase)

rule run_hic_bowtie:
    input: expand("{projectdir}/results/{datatype}/{srr}_bowtie_aligned.bam", projectdir=projectdir, datatype="HiC", srr=hic)

rule run_bisulfite_bowtie:
    input: expand("{projectdir}/results/{datatype}/{srr}_bowtie_aligned.bam", projectdir=projectdir, datatype="bisulfite", srr=bisulfite)

rule run_chipseq_bowtie:
    input: expand("{projectdir}/results/{datatype}/{srr}_bowtie_aligned.bam", projectdir=projectdir, datatype="chipseq", srr=chipseq)


## rules for maping with minimap2
rule run_dnase_minimap2:
    input: expand("{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam", projectdir=projectdir, datatype="dnase", srr=dnase)

rule run_hic_minimap2:
    input: expand("{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam", projectdir=projectdir, datatype="HiC", srr=hic)

rule run_bisulfite_minimap2:
    input: expand("{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam", projectdir=projectdir, datatype="bisulfite", srr=bisulfite)

rule run_chipseq_minimap2:
    input: expand("{projectdir}/results/{datatype}/{srr}_minimap2_aligned.bam", projectdir=projectdir, datatype="chipseq", srr=chipseq)
