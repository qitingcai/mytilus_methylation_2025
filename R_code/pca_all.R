library(vegan)
library(edgeR)
library(tidyverse)
library(mice)
library(dplyr)

setwd("/hb/groups/kelley_lab/tina/mytilus/Rscripts/data/merged_cov/top5")

targets <- read.delim("sample_full_top5_nobaseline.txt", row.names = "sample", stringsAsFactors = FALSE)

meta<-targets %>%
  rownames_to_column(var = "sample_name")
Sample <- row.names(targets)

meta_data<-read.delim("sample_full_top5_nobaseline.txt", row.names = "sample", stringsAsFactors = FALSE)

files <- paste0(Sample,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")
yall <- readBismark2DGE(files, sample.names=Sample)

#pca on original samples

Methylation <- gl(2, 1, ncol(yall), labels = c("Me", "Un"))
Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]

#pca on original samples
TotalLibSize <- yall$samples$lib.size[Methylation=="Me"] +
  +                 yall$samples$lib.size[Methylation=="Un"]
yall$samples$lib.size <- rep(TotalLibSize, each=2)
yall$samples

Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]

prop_meth_matrix <- Me/(Me+Un)

# Ensure column names start with a letter
colnames(prop_meth_matrix) <- make.names(colnames(prop_meth_matrix), unique = TRUE)
colnames(prop_meth_matrix)
imputed_data <- mice(prop_meth_matrix, m = 1, method = 'pmm', seed = 123)

# Complete the data using the first imputed dataset
complete_data <- complete(imputed_data, 1)

write.csv(complete_data, "imp_matrix_merged_40samples.csv")

# Perform PCA using prcomp
imp_pca_result <- prcomp(t(complete_data))

# Create df of pca loadings with imputation
imp_scores <- as.data.frame(imp_pca_result$x)

write.csv(imp_scores, "imp_scores_pca_merged_40samples.csv")

pca_scores<-read.csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/imp_scores_pca_merged_40samples.csv")


check_G_between <- function(x) {
  # Find '.G_'
  if (grepl("\\.G_", x)) {
    return("Gill")
  } else {
    return("Foot")
  }
}

rownames(pca_scores)<-pca_scores[,1]
pca_scores$Tissue <- sapply(row.names(pca_scores), check_G_between)

sample_ids <- rownames(pca_scores)
extracted_letters <-sub("^X(\\w+).*", "\\1", sample_ids)

meta_data<-read.delim("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/sample_full_top5_nobaseline.txt", row.names = "sample", stringsAsFactors = FALSE) %>% 
  na.omit()

# Change the first column to keep only the desired portion
pca_scores[, 1] <- sub("^X(\\w+)\\..*", "\\1", pca_scores[, 1])


pca_data <- data.frame(
  PC1 = pca_scores$PC1,
  PC2 = pca_scores$PC2,
  Label = extracted_letters,
  SampleID = extracted_letters ,
  Tissue = meta_data$tissue# Assuming 'sample_ids' is a vector containing sample IDs
)

# Create lines connecting points with the same letter
line_data <- pca_data %>%
  group_by(Label) %>%
  filter(n() > 1) %>%
  arrange(Label, PC1, PC2)

pca_imputated_merged<- ggplot(pca_data, aes(x = PC1, 
                                              y = PC2, 
                                              color = meta_data$tissue,
                                              shape = meta_data$origin_site)) +
  geom_point(size = 8) +
  geom_line(aes(group = SampleID), color = "black", linetype = 1) +
  #geom_text(aes(label = SampleID), hjust = 1.5, vjust = 1.5, size = 3, check_overlap = TRUE) +
  theme_bw() +
  labs(
    x = "PC1 (5.14%)",
    y = "PC2 (4.71%)",
    color = "Tissue Type",
    shape = "Origin Site"
  )+
  scale_color_manual(
    values = c("F" = "#3bbdcc", "G" = "#e68d4e"),  # Map colors
    labels = c("F" = "Foot", "G" = "Gill")         # Set legend labels
  ) +
  scale_shape_manual(
    values = c("protected" = 16, "exposed" = 15),  # Customize shapes, e.g., 16 = circle, 17 = triangle
    labels = c("protected" ="Protected", "exposed"="Exposed")             # Set shape legend labels
  )+
  theme(
    panel.grid = element_blank()  # Remove grid lines
  )+theme(
    strip.text = element_text(size = 14),         # Increase size of facet strip text
    axis.text.x = element_text(size = 12), # Increase x-axis text size
    axis.text.y = element_text(size = 12),       # Increase y-axis text size
    axis.title = element_text(size = 14),        # Increase axis titles size
    legend.text = element_text(size = 12),       # Increase legend text size
    legend.title = element_text(size = 14),      # Increase legend title size
    panel.grid = element_blank()                 # Remove background grid lines
  )

pca_imputated_merged




