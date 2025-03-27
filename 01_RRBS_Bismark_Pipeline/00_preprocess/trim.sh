#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=trimgalore            # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=9                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=30G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

### Aim of this script: code for paired-end 65bp hard trimming on all samples ###

### Loading trimgalore, -> containing TrimGalore version 0.6.10 ###
module load trimgalore 

### Define the input directory containing the raw FASTQ files ###
INPUT_DIR="/hb/groups/kelley_lab/tina/mytilus/data"

### Define the output directory where trimmed FASTQ files will be saved ###
OUTPUT_DIR="/hb/groups/kelley_lab/tina/mytilus/02_trim/final_trim0603"

### Loop through each pair of FASTQ files ###
for file1 in $INPUT_DIR/*_R1_001.fastq.gz; do
    # Extract sample name and lane number from the file name
    sample=$(basename "$file1" | cut -d'_' -f1-3)
    
    # Loop through each lane (L001 to L004)
    for lane in {001..004}; do
        # Construct filenames for R1 and R2 in the current lane
        current_file1="${sample}_R1_${lane}.fastq.gz"
        current_file2="${sample}_R2_${lane}.fastq.gz"
        
        # Check if the current pair of files exist
        if [ -e "$INPUT_DIR/$current_file1" ] && [ -e "$INPUT_DIR/$current_file2" ]; then
            # Run TrimGalore on the current pair of files
            trim_galore --gzip "$INPUT_DIR/$current_file1" "$INPUT_DIR/$current_file2" --output_dir "$OUTPUT_DIR" --hardtrim5 65 --cores 2 --paired --rrbs --stringency 1 --fastqc_args "--nogroup" #no group show all the basepairs 
        fi
    done
done