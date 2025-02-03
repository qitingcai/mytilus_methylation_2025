#####foot--- run at same time
setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/tissue_specific/foot")
meta_data_foot<-read.delim("foot_metadata.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_foot<- row.names(meta_data_foot)

# List of file names
files_foot <- paste0(Sample_foot,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

# Load and combine the files
combined_coverage <- do.call(rbind, lapply(files_foot[file.exists(files_foot)], function(file) {
  data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  sample_id <- gsub("\\.CpG_report\\.merged_CpG_evidence\\.cov\\.CpG_report\\.merged_CpG_evidence\\.cov$", "", basename(file))
  data$Sample <- sample_id
  return(data)
}))


# Rename columns (example, modify as per your file structure)
colnames(combined_coverage) <- c("Chromosome", "Position", "Strand", "Coverage", "Methylated_Count", "Unmethylated_Count", "Sample")


meth_perc_cov <- combined_coverage %>% 
  group_by(Sample) %>% 
  summarize(Raw_perc = mean(Coverage))

#######calculate methylation percentage from the prop_meth_matrix_gill

prop_meth_perc <- as.data.frame(colMeans(prop_meth_matrix_foot, na.rm = TRUE)*100)
# Add the column names of the original matrix as a new column
prop_meth_perc <- cbind(Sample = rownames(prop_meth_perc), prop_meth_perc)

# Rename the columns for clarity
colnames(prop_meth_perc) <- c("Sample", "Filtered_perc")
prop_meth_perc$Sample <- gsub("-Me", "", prop_meth_perc$Sample)

meta_data_foot<-
   cbind(Sample = rownames(meta_data_foot), meta_data_foot)

combined_meth_perc<-meth_perc_cov %>% 
  left_join(prop_meth_perc, by ="Sample") %>% 
  left_join(meta_data_foot,  by ="Sample")