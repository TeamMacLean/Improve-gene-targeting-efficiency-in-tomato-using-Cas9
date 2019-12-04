#!/bin/bash

## bash script to download SRA files using sra run table from ncbi

sraruntable=$1
outputfolder=$2


for datatype in $(grep -i solanum  $sraruntable  |awk -F "," '{print $2}' | sort -r | uniq); do
	# this loop will look for these datatype : RNA-Seq, OTHER(HiC), DNase-Hypersensitivity, ChIP-Seq, Bisulfite-Seq
	
	for foldername in $(grep -i solanum  $sraruntable | grep $datatype |  sed 's/gs.US,ncbi.public,s3.us-east-1/gs.US-ncbi.public-s3.us-east-1/;s/gs,ncbi,s3/gs-ncbi-s3/;s/ /_/g;s/,/\t/g' | awk 'BEGIN{OFS="--"}{print $1, $2, $17,$19, $24,$28, $30, $31}'  | sed 's/OTHER/HiC/' | grep -v leaf_old | grep -v fruit); do
		#echo creating folder ${outputfolder}/$foldername
		mkdir -p ${outputfolder}/$foldername
		srr=$(echo $foldername | awk -F "--" '{print $1}' )
		cmd="fastq-dump --split-files --gzip --outdir ${outputfolder}/$foldername $srr"
		echo command $cmd
		# run the fastq-dum command

		$cmd

	done
done

