##### Aim of this script: 
##### calculating the average methylation percentages across all CpG sites across samples (for both foot and gill) before and after CpG site filtering #####


############## Foot subset ############## 

setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/tissue_specific/foot")
meta_data_foot<-read.delim("foot_metadata.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_foot<- row.names(meta_data_foot)

### List of file names ###
files_foot <- paste0(Sample_foot,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

### Load and combine the files ###
combined_coverage <- do.call(rbind, lapply(files_foot[file.exists(files_foot)], function(file) {
  data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  sample_id <- gsub("\\.CpG_report\\.merged_CpG_evidence\\.cov\\.CpG_report\\.merged_CpG_evidence\\.cov$", "", basename(file))
  data$Sample <- sample_id
  return(data)
}))

### Rename columns ###
colnames(combined_coverage) <- c("Chromosome", "Position", "Strand", "Coverage", "Methylated_Count", "Unmethylated_Count", "Sample")

### Calculatate raw average methylation proportion for each sample ###
meth_perc_cov <- combined_coverage %>% 
  group_by(Sample) %>% 
  summarize(Raw_perc = mean(Coverage))

####### Calculate filtered average methylation percentage from the prop_meth_matrix_foot (from foot_DM_analyses.R) ####
prop_meth_perc <- as.data.frame(colMeans(prop_meth_matrix_foot, na.rm = TRUE)*100)
# Add the column names of the original matrix as a new column
prop_meth_perc <- cbind(Sample = rownames(prop_meth_perc), prop_meth_perc)

### Rename the columns for clarity ###
colnames(prop_meth_perc) <- c("Sample", "Filtered_perc")
prop_meth_perc$Sample <- gsub("-Me", "", prop_meth_perc$Sample)

### Create dataframw with both raw and filtered average methylation percentage calculations ###
meta_data_foot<-
   cbind(Sample = rownames(meta_data_foot), meta_data_foot)

combined_meth_perc<-meth_perc_cov %>% 
  left_join(prop_meth_perc, by ="Sample") %>% 
  left_join(meta_data_foot,  by ="Sample")


############## Gill subset ############## 

setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/tissue_specific/gill")
meta_data_gill<-read.delim("gill_metadata.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_gill <- row.names(meta_data_gill)

### List of file names ###
files_gill <- paste0(Sample_gill,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

### Load and combine the files ###
combined_coverage <- do.call(rbind, lapply(files_gill[file.exists(files_gill)], function(file) {
  data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  sample_id <- gsub("\\.CpG_report\\.merged_CpG_evidence\\.cov\\.CpG_report\\.merged_CpG_evidence\\.cov$", "", basename(file))
  data$Sample <- sample_id
  return(data)
}))

### Rename columns ###
colnames(combined_coverage) <- c("Chromosome", "Position", "Strand", "Coverage", "Methylated_Count", "Unmethylated_Count", "Sample")

### Calculatate raw average methylation proportion for each sample ###
meth_perc_cov <- combined_coverage %>% 
  group_by(Sample) %>% 
  summarize(Raw_perc = mean(Coverage))


####### Calculate filtered average methylation percentage from the prop_meth_matrix_foot (from gill_DM_analyses.R) ####
prop_meth_perc <- as.data.frame(colMeans(prop_meth_matrix_gill, na.rm = TRUE)*100)
# Add the column names of the original matrix as a new column
prop_meth_perc <- cbind(Sample = rownames(prop_meth_perc), prop_meth_perc)

### Rename the columns for clarity ###
colnames(prop_meth_perc) <- c("Sample", "Filtered_perc")
prop_meth_perc$Sample <- gsub("-Me", "", prop_meth_perc$Sample)

### Create dataframw with both raw and filtered average methylation percentage calculations ###
meta_data_gill<-
   cbind(Sample = rownames(meta_data_gill), meta_data_gill)

combined_meth_perc<-meth_perc_cov %>% 
  left_join(prop_meth_perc, by ="Sample") %>% 
  left_join(meta_data_gill,  by ="Sample")



############## Making box plots of average methylation percentages ############## 
### Need to have a dataframe that combines gill and foot subset information before making the plot ###
meth_file<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/combined_tissue_meth_perc.csv")

### Loading libraries ###
library(ggplot2)
library(tidyr)
library(dplyr)


### data wrangling for plotting ###
reshaped_meth_file <- meth_file  %>%
  pivot_longer(
    cols = c(Raw_perc, Filtered_perc),
    names_to = "Percentage_Type",
    values_to = "Percentage"
  )

Gill<-reshaped_meth_file %>% 
  filter(tissue =="Gill",
         Percentage_Type == "Filtered_perc")

Foot<-reshaped_meth_file %>% 
  filter(tissue =="Foot",
         Percentage_Type == "Filtered_perc")


### Plot average methylation percentage by tissue type ###
ggplot(meth_file, aes(x = tissue, y = Raw_perc, fill = tissue)) + 
  geom_boxplot() +
  labs(
    title = "",
    x = "Tissue Type",
    y = "Global Methylation Level (%)",
    fill = "Tissue Type"
  ) +
  scale_fill_manual( values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e"))+
  theme_classic() +
  theme(
    strip.text = element_text(size = 12),
    panel.grid = element_blank() # Remove background grid lines
  )+theme(
    strip.text = element_text(size = 14),         # Increase size of facet strip text
    axis.text.x = element_text(size = 12), # Increase x-axis text size
    axis.text.y = element_text(size = 12),       # Increase y-axis text size
    axis.title = element_text(size = 14),        # Increase axis titles size
    legend.text = element_text(size = 12),       # Increase legend text size
    legend.title = element_text(size = 14),      # Increase legend title size
    panel.grid = element_blank()                 # Remove background grid lines
  )+geom_jitter(height=0)




### Facet by transplant and origin site ###

meth_file_long <- meth_file %>%
  pivot_longer(cols = c(origin_site, transplant_site),  
               names_to = "site_type",                   
               values_to = "site") 

ggplot(meth_file_long, aes(x = site, y = Raw_perc, fill = tissue)) + 
  geom_boxplot() +
  facet_wrap(~ site_type, scales = "free_x") +  # Facet by site_type (either origin or transplant)
  labs(
    title = " ",
    x = " ",
    y = "Global Methylation Level (%) ",
    fill = "Tissue Type"
  ) +
  scale_fill_manual(values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e")) +  # Color mapping for tissue
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12),    # Increase facet title size
    panel.grid = element_blank()  # Remove background grid lines
  )+theme(
    strip.text = element_text(size = 14),         # Increase size of facet strip text
    axis.text.x = element_text(size = 12), # Increase x-axis text size
    axis.text.y = element_text(size = 12),       # Increase y-axis text size
    axis.title = element_text(size = 14),        # Increase axis titles size
    legend.text = element_text(size = 12),       # Increase legend text size
    legend.title = element_text(size = 14),      # Increase legend title size
    panel.grid = element_blank()                 # Remove background grid lines
  )+geom_jitter(height=0)



