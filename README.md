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
│   │   │      ├── concat_lanes.sh
│   │   │      ├── concat_fastqc.sh
│   │   │      └── concat_multiqc.sh
│   │   ├── Bismark_01
│   │   │      ├── genome_download.sh
│   │   │      ├── genome_preparation.sh
│   │   │      ├── alignment.sh
│   │   │      ├── methy_extractor.sh
│   │   └──    └── merge_strand.sh
│   │  
│   ├── 02_R_code
│   │   ├── foot_DM_analysis.R
│   │   ├── gill_DM_analysis.R
│   │   ├── pca_all.R
│   │   └── average_meth.R
│   └── Readme.md
│
└── Data
    ├── Bismark coverage files (n=40)



```
