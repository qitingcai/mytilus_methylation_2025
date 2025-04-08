library(vegan)
library(edgeR)
library(tidyverse)
library(mice)
library(dplyr)

setwd("/hb/groups/kelley_lab/tina/mytilus/Rscripts/data/merged_cov/top5")

meta_data <- read.delim("sample_full_top5_nobaseline.txt", row.names = "sample", stringsAsFactors = FALSE)

meta<-meta_data< %>%
  rownames_to_column(var = "sample_name")
Sample <- row.names(meta_data<)

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

#write.csv(complete_data, "imp_matrix_merged_40samples.csv")

# Perform PCA using prcomp
imp_pca_result <- prcomp(t(complete_data))

# Create df of pca loadings with imputation
imp_scores <- as.data.frame(imp_pca_result$x)

#write.csv(imp_scores, "imp_scores_pca_merged_40samples.csv")

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

meta_data<-read.delim("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/sample_full_top5_nobaseline_pca.txt", row.names = "sample", stringsAsFactors = FALSE) %>% 
  na.omit()

# Change the first column to keep only the desired portion
pca_scores[, 1] <- sub("^X(\\w+)\\..*", "\\1", pca_scores[, 1])

### Merging metadata and pca scores data ###
pca_data <-meta_data %>% 
  left_join(pca_scores,
            by=c("meth"="X",
                 "tissue"="Tissue"))

### Create lines connecting points with the same letter ###
line_data <- pca_data %>%
  group_by(meth) %>%
  filter(dplyr::n() > 1) %>%
  arrange(meth, PC1, PC2)

pca_imputated_merged <- ggplot(pca_data, aes(
  x = PC1,
  y = PC2,
  shape = GROUP_CODE,
  fill = tissue
)) +
  geom_point(size = 10, alpha = 0.85, stroke = 0.2, color = "black") +
  geom_line(aes(group = meth), color = "black", linetype = 1) +
  scale_shape_manual(
    values = c(
      "protected_protected" = 21,
      "protected_exposed" = 22,
      "exposed_exposed" = 23,
      "exposed_protected" = 24
    ),
    labels = c(
      "protected_protected" = "Protected_Protected",
      "protected_exposed" = "Protected_Exposed",
      "exposed_exposed" = "Exposed_Exposed",
      "exposed_protected" = "Exposed_Protected"
    )
  ) +
  scale_fill_manual(
    values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend()
  ) +
  labs(
    x = "PC1 (5.14%)",
    y = "PC2 (4.71%)",
    fill = "Tissue Type",
    shape = "Treatment"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pca_imputated_merged


### Assessing relationship between PC1 and final shell length ###

### Add length data ###
length_data<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/sample_librarysize_length_top5.csv")

### Run a regression of PC1 and final shell length ###

#reformatting
pca_scores$sample_id<-rownames(pca_scores) %>% 
  gsub("^X", "", .) %>%        # Remove the initial "X"
  gsub("_.*$", "", .) %>% 
  gsub("\\.", "-", .)  

colnames(length_data)[1] <- "sample_name"

### combining PCA scores and sample information (length data) ###
PC <- length_data %>% 
  left_join(pca_scores, by = c("sample ID"="sample_id"))

PC$`sample ID` <- sub("-.*", "", PC$`sample ID`)

pc1_length<-ggplot(PC, aes(x = `LENGTH_FINAL (mm)`, y = PC1)) + #plotting shell length and PC1 loadings
  geom_point(aes(color = Tissue.x),
             size=10) +  # Keep points colored by tissue
  geom_smooth(method = "lm", size = 0.8, alpha = 0.1, color = "black") +  # Single regression line for all shell length and PC1
  geom_line(aes(group = `sample ID`), color = "black", linetype = 1) + 
  labs(
    y = "PC1", #label y axis
    x = "Final shell length (mm)", #label x axis
    color = "Tissue Type" #specify color of points
  ) + scale_color_manual(values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e"),
                         labels = c("F" = "Foot", "G" = "Gill"))+ #specify color of points
  theme_bw() +
  theme( strip.text = element_text(size = 20),  # Facet strip text size
         axis.text.x = element_text(size = 20), # X-axis text size
         axis.text.y = element_text(size = 20), # Y-axis text size
         axis.title = element_text(size = 20),  # Axis title size
         legend.text = element_text(size = 20), # Legend text size
         legend.title = element_text(size = 20),# Legend title size
         panel.grid = element_blank()           # Remove grid lines
  )

pc1_length

library(lmerTest)
library(lme4)

### Running GLMM to see if length correlates with PC1 loading ###
lm_length<-lmer(PC1~`LENGTH_FINAL (mm)` + Tissue.x + (1|`sample ID`), data=PC)
summary(lm_length)





