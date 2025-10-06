# mytilus_methylation_2025
Code repository for mussel DNA methylation project

**01_RRBS_Bismark_Pipeline**: this folder includes scripts used to preprocess raw RRBS reads and generate CpG methylation calls via the Bismark pipeline
1. Preprocess_00 (quality assessment, trimming, concatenating lanes, intron annotation on gff file)
2. Bismark_01 (downloading & preparing reference genome, alignment, methylation extraction, merging strand information)

**02_Rcode**: this folder includes scripts used to run differential methylation analyses in R. Analyses were ran separately for foot and gill sample subsets. Also include scripts to make plots and run GO analyses.


**Data**: contains all coverage files used for DM analyses, genome annotatation file (with introns annotated), GO terms

**Outputs**: contains analyses outputs

### File structure
```
├── Mytilus_methylation_2025
│   ├── 01_RRBS_Bismark_Pipeline
│   │   ├── Preprocess_00
│   │   │      ├── pretrim_fastqc_00.sh
│   │   │      ├── pretrim_multiqc_01.sh
│   │   │      ├── trim_02.sh
│   │   │      ├── concat_lanes_03.sh
│   │   │      ├── concat_fastqc_04.sh
│   │   │      ├── concat_multiqc_05.sh
│   │   │      └── add_intron.sh
│   │   ├── Bismark_01
│   │   │      ├── genome_download_00.sh
│   │   │      ├── genome_preparation_01.sh
│   │   │      ├── alignment_02.sh
│   │   │      ├── methy_extractor_03.sh
│   │   └──    └── merge_strand_04.sh
│   │  
│   ├── 02_R_code
│   │   ├── Figures.rmd
│   │   ├── Fishers_tests.rmd
│   │   ├── GOMWU.R
│   │   ├── average_lfc.R
│   │   ├── final_foot_DM_analysis.rmd
│   │   ├── final_gill_DM_analysis.rmd
│   │   ├── interaction_model_foot.rmd
│   │   ├── interaction_model_gill.rmd
│   │   └── model_correlation.rmd
│   ├── Data 
│   │   ├── foot_coverage_files (n=20 + foot_metadata.txt)
│   │   ├── gill_coverage_files (n=20 + gill_metadata.txt)
│   │   ├── GCF_021869535.1_xbMytCali1.0.p_gene_ontology.gaf
│   │   ├── new_genomic_intron.gff
│   │   ├── R input files (.csv)
│   ├── Outputs
|   │   ├── Mbias_output (mbias_report.html)
|   │   ├── multiqc_outputs (trimmed + pretrimmed reports)
|   │   ├── imp_scores_pca_merged_40samples.csv
│   └── └── bismark_output.csv
│   
└── README.MD


```
