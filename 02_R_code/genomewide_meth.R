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


### test whether genome wide DNA methylation proportion differ by tissue type or treatment group, plot the dara ###
### loading data ###
meth_prop_all<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/percentage_allsample_genomewide.csv")
meth_prop_all$`Pre-filtered mean global methylation (%)`<-meth_prop_all$`Pre-filtered mean global methylation (%)`/100
meth_prop_all$`Post-filtered mean global methylation (%)`<-meth_prop_all$`Post-filtered mean global methylation (%)`/100

### Run GLMM ###
library(glmmTMB)
glm<-glmmTMB(`Pre-filtered mean global methylation (%)`~Tissue+`Origin site` + `Transplant site` + (1|Sample),  data=meth_prop_all, family=beta_family(link="logit"))
summary(glm)

glm_post<-glmmTMB(`Post-filtered mean global methylation (%)`~Tissue+`Origin site` + `Transplant site` + (1|Sample),  data=meth_prop_all, family=beta_family(link="logit"))
summary(glm_post)


### jitter with mean and SD ###
p <- ggplot(meth_prop_all_pre, aes(x = `Transplant site`, y = `Pre-filtered mean global methylation (%)`, fill = Tissue)) + 
  # Jittered points with custom appearance
  geom_jitter( size = 8, alpha = 0.4, 
               stroke = 0, shape = 21, 
               position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.6)) + 
  # Adding mean and standard deviation as points with error bars
  stat_summary(
    aes(color = Tissue),  # Correctly map color to Tissue
    fun.data = "mean_sdl", 
    fun.args = list(mult = 1), 
    geom = "pointrange", 
    position = position_dodge(width = 0.6), 
    size = 2,
    shape=18,
    alpha=0.9,
    show.legend = FALSE
  ) + 
  facet_wrap(~`Origin site`) +  # Facet by site_type
  theme(strip.text = element_text(size = 20, face = "bold"),  
        strip.background = element_blank(),  # Removes the gray background
        strip.placement = "outside")+
  scale_fill_manual(values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e")) +  # Custom fill colors for Tissue
  scale_color_manual(values = c("Foot" = "#2c8e99", "Gill" = "#d4562a")) +  # Custom point color for Tissue
  labs(y = "Methylation Proportion") +  # Axis labels and title
  #scale_y_continuous(limits = c(0.1, 0.3)) + 
  # Set y-axis limits from 0 to 1
  theme_few() +  # Minimal theme
  theme(
    strip.text = element_text(size = 20),  # Facet label size
    axis.text = element_text(size = 20),   # Axis text size
    axis.title = element_text(size = 20),  # Axis title size
    legend.title = element_blank(),
   legend.text = element_text(size = 20))

# Print the plot
print(p)

