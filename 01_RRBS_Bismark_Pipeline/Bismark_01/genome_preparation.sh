#!/bin/bash

#SBATCH --partition=256x44               # Partition/queue to run on
#SBATCH --time=2-00:00:00                # Max time for job to run
#SBATCH --job-name=genome_prep                  # Name for job (shows when running squeue)
#SBATCH --mail-type=ALL                  # Mail events(NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=qcai17@ucsc.edu    # Where to send mail
#SBATCH --ntasks=1                       # Number of tasks to run
#SBATCH --cpus-per-task=10                # Number of CPU cores to use per task
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --mem=40G                        # Ammount of RAM to allocate for the task
#SBATCH --output=slurm_%j.out            # Standard output and error log
#SBATCH --error=slurm_%j.err             # Standard output and error log
#SBATCH --no-requeue                     # don't requeue the job upon NODE_FAIL

#### Aim of this script: 
### Bisulfite convert reference genome Mytilus Californianus for downstream Bismark analyses ###

######################### Setting up the environment #########################
### Code for downloading Bismark to home directory, --> containing version v0.24.2 ###
# wget https://github.com/FelixKrueger/Bismark/archive/master.zip 
# unzip master.zip

### Code for downloading Bismark as a Conda environment, --> containing version v0.24.2 ###
# module load miniconda3 
# conda create -n bismark 
# conda install bioconda::bismark 
# conda activate bismark

### Code for downloading Bowtie2 to home directory, --> containing version 2.5.4 ###
# wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.4/bowtie2-2.5.4-sra-linux-x86_64.zip
# unzip download

################################################################################

# Set path to find bismark, bowtie2 and folder that contains the reference genome
bismark_dir="/hb/groups/kelley_lab/tina/mytilus/Bismark-master"
bowtie2_dir="/hb/groups/kelley_lab/tina/mytilus/bowtie2/bowtie2-2.5.4-sra-linux-x86_64/"
genome_folder="/hb/groups/kelley_lab/tina/mytilus/ref_genome/GCF_021869535.1/"

### Running Bismark genome preparation step ###
### Bismark will create two individual folders within this directory, one for a C->T converted genome and the other one for the G->A converted genome. ###
${bismark_dir}/bismark_genome_preparation --verbose --parallel 28 --path_to_aligner ${bowtie2_dir} ${genome_folder}