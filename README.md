# mytilus_methylation_2025
Code repository for mussel DNA methylation project

**01_RRBS_Bismark_Pipeline**: this folder includes scripts used to preprocess raw RRBS reads and generate CpG methylation calls via the Bismark pipeline
1. Preprocess_00 (quality assessment, trimming, concatenating lanes)
2. Bismark_01 (downloading & preparing reference genome, alignment, methylation extraction, merging strand information)

**02_Rcode**: this folder includes scripts used to run differential methylation analyses in R. Analyses were ran separately for foot and gill sample subsets





### File structure
```
├── MYTILUS_METHYLATION_2025
│   ├── 01_RRBS_Bismark_Pipeline
│   │   ├── Preprocess_00
│   │   │      ├── pretrim_fastqc_00.sh
│   │   │      ├── pretrim_multiqc_01.sh
│   │   │      ├── trim_02.sh
│   │   │      ├── concat_lanes_03.sh
│   │   │      ├── concat_fastqc_04.sh
│   │   │      └── concat_multiqc_05.sh
│   │   ├── Bismark_01
│   │   │      ├── genome_download_00.sh
│   │   │      ├── genome_preparation_01.sh
│   │   │      ├── alignment_02.sh
│   │   │      ├── methy_extractor_03.sh
│   │   └──    └── merge_strand_04.sh
│   │  
│   ├── 02_R_code
│   │   ├── foot_DM_analysis.R
│   │   ├── gill_DM_analysis.R
│   │   ├── pca_all.R
│   │   └── average_meth.R
│   └── Readme.md
│
└── Data
    ├── foot_coverage_files (n=20 + foot_metadata.txt)
    └── gill_coverage_files (n=20 + gill_metadata.txt)


```
