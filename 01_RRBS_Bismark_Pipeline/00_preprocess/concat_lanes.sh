#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=3-00:00:00                # Max time for job to run
#SBATCH --job-name=concat                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu      # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20               # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Amount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # Don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-700]                  # Array job, adjust the range as needed

### Aim of this script: 
### concatenating lane 01 to 04 for each sample for both forward and reverse reads 

### Load the sample name based on second column in sample.txt (metadata file) ###
LINE=$(sed -n "${SAMPLE_ID}p" /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/sample.txt)
sample=$(echo ${LINE} | awk '{ print $2; }')

### Define the input and output file paths ###
input_path="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603"
output_path="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603/concat"

### create output directory ###
mkdir -p $output_path

### Concatenate R1 files and save the new merged files ###
cat ${input_path}/${sample}_L001_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L002_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L003_R1_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L004_R1_001.65bp_5prime.fq.gz > ${output_path}/${sample}_R1_merged.fq.gz

### Concatenate R2 files and save the new merged files ###
cat ${input_path}/${sample}_L001_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L002_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L003_R2_001.65bp_5prime.fq.gz \
    ${input_path}/${sample}_L004_R2_001.65bp_5prime.fq.gz > ${output_path}/${sample}_R2_merged.fq.gz


echo "Finished processing sample: $sample"