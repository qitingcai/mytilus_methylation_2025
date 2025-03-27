#!/bin/bash

#SBATCH --partition=128x24               # Partition/queue to run on
#SBATCH --time=1-00:00:00                # Max time for job to run
#SBATCH --job-name=multiqc                # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=8                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=10G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

### Set directory which contains fastqc result files ###
cd /hb/groups/kelley_lab/tina/mytilus/01_fastqc

### Loading multiqc, -> contains multiqc, version 1.25.2 ###
module load miniconda3
conda activate multiqc 

### Running multiqc on all fastqc files and create aggregated results (.html) for easier viewing ###
multiqc -v -o multiqc/ /hb/groups/kelley_lab/tina/mytilus/01_fastqc/results/