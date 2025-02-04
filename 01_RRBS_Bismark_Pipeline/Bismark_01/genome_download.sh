
#### Aim of this script: 
### Downloadig reference genome Mytilus Californianus https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_021869535.1/ ###

### create a conda environment and download the NCBI Datasets command-line tools (CLI) to access sequence data ###
### loading the CLI ###
module load  miniconda3
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

#### Code to download the NCBI reference genome ###
datasets download genome accession GCF_021869535.1 --include gff3,rna,cds,protein,genome,seq-report

#### need to change sequence extension from .fna to .fa/.fasta to run Bismark downstream ###
cd GCF_021869535.1/
mv GCF_021869535.1_xbMytCali1.0.p_genomic.fna  GCF_021869535.1_xbMytCali1.0.p_genomic.fa