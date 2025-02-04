#!/bin/bash

#SBATCH --partition=Build               # Partition/queue to run on
#SBATCH --time=0-05:00:00                # Max time for job to run
#SBATCH --job-name=methylation_extract   # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu      # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=20               # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=40G                        # Amount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

export PATH="$PATH:/hb/groups/kelley_lab/tina/mytilus/04_methylation_extractor/samtools-1.9"

bismark_dir="/hb/groups/kelley_lab/tina/mytilus/Bismark-master"
map_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/F"

# replace the path with Gill sample subset
# map_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map/G"

# Load necessary modules
module load bismark

# Define the BAM files to process
bam_files=$(ls ${map_dir}/*/*.bam)

### Run Bismark methylation extractor ###

#the output of the methylation extractor can be transformed into a bedGraph and coverage file using the option --bedGraph (see also --counts)
# Context-dependent methylation output files (--comprehensive option) This will default to three output files (CpG-context, CHG-context and CHH-context)
# for paired-end reads --no_overlap (on by default): specifying this option will extract the methylation calls of overlapping parts in the middle of paired-end reads only once (using the calls from the first read which is presumably the one with a lowest error rate).
#-buffer_size <string> This allows you to specify the main memory sort buffer when sorting the methylation information

${bismark_dir}/bismark_methylation_extractor --bedGraph --counts -p --comprehensive --no_overlap --multicore 28 --buffer_size 75% \
--gzip -o /hb/groups/kelley_lab/tina/mytilus/04_methylation_extractor/final_minscore0.6 ${bam_files}