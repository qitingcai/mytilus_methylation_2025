#!/bin/bash

#SBATCH --partition=128x24                # Partition/queue to run on
#SBATCH --time=0-12:00:00                # Max time for job to run
#SBATCH --job-name=bismark                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=4                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=40G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL
#SBATCH --array=[1-200]                 # array job

### Aim of this script: mapping all the preprocessed (trimmed + concatenated) RRBS reads to the bisulfite converted reference genome ####

### Loading Bismark, --> contain  Bismark Version: v0.24.2 ###
module load parallel miniconda3
conda activate bismark

### Load meta file that contains sample name, tissue type, lane (all just L001 because they were concatenated), and strand (R1 or R2) ###
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/sample_final.txt)
sample=$(echo ${LINE} | awk '{ print $2; }')
tissue=$(echo ${LINE} | awk '{ print $3; }')
lane=$(echo ${LINE} | awk '{ print $4; }')
strand=$(echo ${LINE} | awk '{ print $5; }')

echo "running Bismark for sample: ${sample} (${tissue}, ${lane})"

### Defining path to Bowtie2, reference genome folder, and directory that contained all the concatenated, trimmed RRBS reads ###
genome_folder="/hb/groups/kelley_lab/tina/mytilus/ref_genome/GCF_021869535.1/" # this directoy contains both the unmodified genome (as .fa or .fasta files) as well as the two bisulfite genome subdirectories 
trimmed_dir="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603/concat" # this directory contains the trimmed, concatenated reads for all samples
bismark_dir="/hb/groups/kelley_lab/tina/mytilus/Bismark-master" # where bismark is downloaded

### Create output directory based on tissue and sample ###
mkdir -p final_map/${tissue}/${sample}
cd final_map/${tissue}/${sample}

### Running paired-end alignment via Bismark ###
# 1. relax default stringent setting from (L,0,-0.2)
# 2. default is directional BS-seq libraries
# 3. specify paired-end alignment

${bismark_dir}/bismark -p 4 ${genome_folder} \
--gzip -score_min L,0,-0.6 \
-1 ${trimmed_dir}/${sample}_R1_merged.fq.gz \
-2 ${trimmed_dir}/${sample}_R2_merged.fq.gz
