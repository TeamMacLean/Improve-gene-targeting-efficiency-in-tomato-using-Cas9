#!/usr/bin/env snakemake
import os
import glob

configfile: "config.yaml"

reference = config["REFERENCE"]
index = config["REFERENCE"]
datafolder = config["datafolder"]
projectdir = config["projectdir"]

# get files from directories in datafolder path
directories_in_datafolder = [directory for directory in glob.glob(datafolder + "/*")]

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

print(datafolder_structure)
# get the folders and files in a dictionary datatype
#
# directory, SRA = glob_wildcards(datafolder + "/{dir}/{srr}_1.fastq.gz")
#
# for x in zip(directory, SRA):
#     print(x)

# rule rnaseq_align_paired:
# 	"""
# 	trinity align reads to reference
# 	"""
# 	input:
# 		R1=lambda wildcards: datafolder_structure[wildcards.srr]["R1"],
#         R2=lambda wildcards: datafolder_structure[wildcards.srr]["R2"],
# 	output: "{projectdir}/results/rnaseq/{srr}_paired_aligned.bam"
# 	shell: "echo tophat2 {input.R1} {input.R2}"

rule bowtie2_build_index:
    input: config['REFERENCE']
    message: "Building index for the reference sequence"
    threads: 4
    shell: "bowtie2-build -f {input} {input}"

rule bowtie2_paired_align:
    input:
        R1=lambda wildcards: datafolder_structure[wildcards.srr]["R1"],
        R2=lambda wildcards: datafolder_structure[wildcards.srr]["R2"]
    output:
        temp("{projectdir}/results/{srr}_paired_bowtie_aligned.sam")
    shell:
        "bowtie2 --no-unal --no-discordant --rgid {wildcards.srr} -x {index} -1 {input.R1} -2 {input.R2} -S {output}"

rule bowtie2_unpaired_align:
    input:"{dir}/{srr}.fastq.gz"
    output: temp("{projectdir}/{results}/{srr}_bowtie_align.sam")
    message: "Aligning single end data with bowtie2 "
    threads: 2
    shell: "bowtie2 --rg-id {wildcards.srr} -x {index} -U {input} -S {output}"

rule run_rnaseq:
    input: expand("{projectdir}/results/{srr}_paired_bowtie_aligned.sam", projectdir=projectdir, srr=datafolder_structure.keys() )
