## Improve-gene-targeting-efficiency-in-tomato-using-Cas9

Improve gene targeting efficiency in tomato using Cas9

## Requirements

1) snakemake-5.5.3
2) python-3.6.1
3) samtools-1/9
4) tophat-2.1.1
5) bowtie2-2.3.5


## How do I run

Run rnaseq data alignment 

```
sbatch --mem 100 -o rnaseqmaster.log -J rnaseq --wrap "snakemake -s scripts/analysis.smk  --jobs 4 -p --verbose --latency-time 60  run_rnaseq --cluster-config cluster.json --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  "
```

Run dnase seq data alignment

```
sbatch --mem 100 -o dnasemaster.log -J dnase  --wrap "snakemake -s scripts/analysis.smk  --jobs 4 -p --verbose --latency-time 60  run_dnase --cluster-config cluster.json --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  "
```

run bisulfite seq data alignment

```
sbatch --mem 100 -o bisulfitemaster.log -J bisulfite --wrap "snakemake -s scripts/analysis.smk  --jobs 4 -p --verbose --latency-time 60 bisulfite --cluster-config cluster.json --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  "
```

run chipseq data alignment

```
sbatch --mem 100 -o chipseqmaster.log -J chipseq --wrap "snakemake -s scripts/analysis.smk  --jobs 4 -p --verbose --latency-time 60  run_chipseq --cluster-config cluster.json --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  "
```
