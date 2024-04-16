# bulk-rnaseq

This pipeline is largely based on the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline. The main purpose of this repository is to provide a tutorial for running the nf-core/rnaseq pipeline on the MSU HPCC using SLURM as the job executor.

This is a workflow to quantify transcript-level RNASeq abundance, detect differentially expressed genes (DEGs), and perform gene set enrichment analysis ('GSEA'). This pipeline uses `STAR` to map FASTQ reads to the reference genome, `Salmon` to perform BAM-level quantification, and `DESeq2` to normalize counts and detect DEGs.

## Resource Requirements
STAR typically uses around 38GB of RAM. It is recommended to run this workflow on an HPC.

## Steps
1. Create a directory for your your analysis (ex: /mnt/home/$username/rnaseq).
2. Make a samplesheet table in csv format. (ex: samplesheet.csv)
3. Make a nextflow configuration file. (ex: nextflow.config)
4. Download the reference genome and gtf files into your directory.
5. Write a bash script to run the pipeline using SLURM. (ex: run_rnaseq.sb)
6. Run the pipeline from your rnaseq directory. (ex: sbatch run_rnaseq.sb)

## Create a directory for your rnaseq analysis
Login to your HPCC account using OnDemand. Navigate to your home directory by clicking 'Files' in the navigation bar. Select 'Home Directory'.

Create a directory for your analysis by clicking 'New Directory'. Name your directory (ex: rnaseq). Navigate to the newly created rnaseq directory.

## Make a samplesheet table
In your rnaseq directory, click 'New File'. Name the file 'samplesheet.csv'. Click the `⋮` symbol and select edit. Create the samplesheet table, for example:
```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,/path/to/fastq/files/CONTROL_REP1_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP1_R2_001.fastq.gz,auto
CONTROL_REP2,/path/to/fastq/files/CONTROL_REP2_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP2_R2_001.fastq.gz,auto
CONTROL_REP3,/path/to/fastq/files/CONTROL_REP3_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP3_R2_001.fastq.gz,auto
```
Save the samplesheet.csv file and return to your rnaseq directory.

## Make a Nextflow configuration file to use SLURM as the process executor
In your rnaseq directory, click `New File`. Name the file 'nextflow.config'. Click the `⋮` symbol and select edit. Create the Nextflow config file:
```
process {
    executor = 'slurm'
}
```
Save the nextflow.config file and return to your rnaseq directory.

## Download the reference genome and gtf files in your rnaseq directory
In your rnaseq directory, click `Open in Terminal` to enter a development node. Download the most recent genome primary assembly and gtf for the organism of interest from [Ensembl](https://ensembl.org/). This may take more than a few minutes. For example:
```
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
```
Note: If you do not download the genome and gtf for the organism that the RNA-seq data is derived from, this pipeline will not work correctly.

## Write a bash script to run the pipeline using SLURM
In your rnaseq directory, click `New File`. Name the file 'run_rnaseq.sb;. Write the bash script, using #SBATCH directives to set resources, for example:
```
#!/bin/bash

#SBATCH --job-name=$jobname_rnaseq
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --cpus-per-task=8

cd /mnt/home/$username/rnaseq
module load Nextflow/23.10.0

nextflow pull nf-core/rnaseq
nextflow run nf-core/rnaseq -r 3.14.0 --input ./samplesheet.csv  -profile singularity --outdir ./rnaseq_results --fasta ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf ./Homo_sapiens.GRCh38.108.gtf.gz -work-dir /mnt/scratch/$username/rnaseq_work -c ./nextflow.config
```
Save the run_rnaseq.sb file and return to your rnaseq directory.

## Run the pipeline
In your rnaseq directory, click `Open in Terminal` to enter a development node. Run jobs on the SLURM cluster:
```
sbatch run_rnaseq.sb
```
Check the status of your pipeline job:
```
squeue -u $username
```
Note: replace $username with your username
