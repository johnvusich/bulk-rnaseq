# Bulk RNA-seq

This pipeline is largely based on the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline. The main purpose of this repository is to provide a tutorial for running the nf-core/rnaseq pipeline on the MSU HPCC using SLURM as the job executor.

This is a workflow to quantify transcript-level RNASeq abundance, detect differentially expressed genes (DEGs), and perform gene set enrichment analysis ('GSEA'). This pipeline uses `STAR` to map FASTQ reads to the reference genome, `Salmon` to perform BAM-level quantification, and `DESeq2` to normalize counts and detect DEGs.

## Resource Requirements
STAR typically uses around 38GB of RAM. It is recommended to run this workflow on an HPC.

## Steps
1. Create a directory for your your analysis (ex: $HOME/rnaseq).
2. Make a samplesheet table for pre-processing in csv format. (ex: samplesheet.csv)
3. Make a nextflow configuration file. (ex: nextflow.config)
4. Download the reference genome and gtf files into your directory.
5. Write a bash script to run the pre-processing pipeline using SLURM. (ex: run_rnaseq.sb)
6. Run the pre-processing pipeline from your rnaseq directory. (ex: sbatch run_rnaseq.sb)
7. Make a samplesheet table for differential expression analysis in csv format. (ex: DE_samplesheet.csv)
8. Write a bash script to run the differential expression pipeline using SLURM. (ex: run_differential.sb)

## Create a directory for your rnaseq analysis
Login to your HPCC account using OnDemand. Navigate to your home directory by clicking 'Files' in the navigation bar. Select 'Home Directory'.

Create a directory for your analysis by clicking 'New Directory'. Name your directory (ex: rnaseq). Navigate to the newly created rnaseq directory.

Make sure to upload your data (FASTQ files) into this directory.

## Make a samplesheet table for RNA-seq pre-processing
In your rnaseq directory, click 'New File'. Name the file 'samplesheet.csv'. Click the `⋮` symbol and select edit. Create the samplesheet table FOR YOUR DATA. Below is a template:
```
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,/path/to/fastq/files/CONTROL_REP1_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP1_R2_001.fastq.gz,auto
CONTROL_REP2,/path/to/fastq/files/CONTROL_REP2_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP2_R2_001.fastq.gz,auto
CONTROL_REP3,/path/to/fastq/files/CONTROL_REP3_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP3_R2_001.fastq.gz,auto
TREATED_REP1,/path/to/fastq/files/TREATED_REP1_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP1_R2_001.fastq.gz,auto
TREATED_REP2,/path/to/fastq/files/TREATED_REP2_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP2_R2_001.fastq.gz,auto
TREATED_REP3,/path/to/fastq/files/TREATED_REP3_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP3_R2_001.fastq.gz,auto
```
Note: the samplesheet above is an example to show the required format. You will still need to make a samplesheet FOR YOUR DATA.

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

cd $HOME/rnaseq
module load Nextflow/23.10.0

nextflow pull nf-core/rnaseq
nextflow run nf-core/rnaseq -r 3.14.0 --input ./samplesheet.csv  -profile singularity --outdir ./rnaseq_results --fasta ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf ./Homo_sapiens.GRCh38.108.gtf.gz -work-dir $SCRATCH/rnaseq_work -c ./nextflow.config
```
Note: above is an example script. You will need to write one specific to YOUR DATA using this template.

Save the run_rnaseq.sb file and return to your rnaseq directory.

## Run the RNA-seq pipeline
In your rnaseq directory, click `Open in Terminal` to enter a development node. Run jobs on the SLURM cluster:
```
sbatch run_rnaseq.sb
```
Check the status of your pipeline job:
```
squeue -u $username
```
Note: replace $username with your username

# Differential Expression Analysis

## Make a samplesheet table for differential expression analysis
In your rnaseq directory, click 'New File'. Name the file 'DE_samplesheet.csv'. Click the `⋮` symbol and select edit. Create the samplesheet table FOR YOUR DATA. Below is a template:
```
sample,fastq_1,fastq_2,condition,replicate,batch
CONTROL_REP1,/path/to/fastq/files/CONTROL_REP1_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP1_R2_001.fastq.gz,control,1,
CONTROL_REP2,/path/to/fastq/files/CONTROL_REP2_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP2_R2_001.fastq.gz,control,2,
CONTROL_REP3,/path/to/fastq/files/CONTROL_REP3_R1_001.fastq.gz,/path/to/fastq/files/CONTROL_REP3_R2_001.fastq.gz,control,3,
TREATED_REP1,/path/to/fastq/files/TREATED_REP1_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP1_R2_001.fastq.gz,treated,1,
TREATED_REP2,/path/to/fastq/files/TREATED_REP2_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP2_R2_001.fastq.gz,treated,2,
TREATED_REP3,/path/to/fastq/files/TREATED_REP3_R1_001.fastq.gz,/path/to/fastq/files/TREATED_REP3_R2_001.fastq.gz,treated,3,
```
Note: the samplesheet above is an example to show the required format. You will still need to make a samplesheet FOR YOUR DATA.

Save the samplesheet.csv file and return to your rnaseq directory.

## Make a contrasts file for differential expression analysis
In your rnaseq directory, click 'New File'. Name the file 'contrasts.csv'. Click the `⋮` symbol and select edit. Create the contrasts file FOR YOUR DATA. Below is a template:
```
id,variable,reference,target,blocking
condition_control_treated,condition,control,treated,
```
Note: the contrasts file above is an example to show the required format. You will still need to make a contrasts file FOR YOUR DATA.

Save the contrasts.csv file and return to your rnaseq directory.

## Write a bash script to run the differential expression analysis using SLURM
In your rnaseq directory, click `New File`. Name the file 'run_differential.sb;. Write the bash script, using #SBATCH directives to set resources, for example:
```
#!/bin/bash

#SBATCH --job-name=$jobname_differential
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --cpus-per-task=8

cd $HOME/rnaseq
module load Nextflow/23.10.0

nextflow pull nf-core/differentialabundance
nextflow run nf-core/differentialabundance -r 1.5.0 --input ./samplesheet.csv --contrasts ./contrasts.csv --matrix ./rnaseq_results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv --gtf ./Homo_sapiens.GRCh38.108.gtf.gz --outdir ./differential_results -profile singularity -work-dir $SCRATCH/differential_work -c ./nextflow.config
```
Note: above is an example script. You will need to write one specific to YOUR DATA using this template.

Save the run_rnaseq.sb file and return to your rnaseq directory.

## Run the differential expression analysis pipeline
In your rnaseq directory, click `Open in Terminal` to enter a development node. Run jobs on the SLURM cluster:
```
sbatch run_differential.sb
```
Check the status of your pipeline job:
```
squeue -u $username
```
Note: replace $username with your username
