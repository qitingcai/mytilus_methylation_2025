##### Script to create Figure 2 #####
### load libraries
library(edgeR)
library(mice)
library(tidyverse)
library(systemfonts)
library(Biobase)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(fields)
library(lmerTest)
library(patchwork)

###### Imputation - run on Hummingbird #####

setwd("/hb/groups/kelley_lab/tina/mytilus/Rscripts/data/merged_cov/outlier_detection")

### Load metadata file that include sample names, load all coverage files into R ###
meta_data_foot<-read.delim("allsamples.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_foot <- row.names(meta_data_foot)
files_foot <- paste0(Sample_foot,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

meta_data_foot_aqm <- meta_data_foot
row.names(meta_data_foot_aqm) <- gsub("-", ".", paste0("X", row.names(meta_data_foot_aqm),".Me"))

### EdgeR function readBismark2DGE reads all the files and collates the counts for all the sample into one data object ###
yall <- readBismark2DGE(files_foot, sample.names=Sample_foot)

### Check dimension and the count matrix ###
dim(yall)
head(yall$counts)

### Save a data frame for downstream analyses ###
yall_df<-as.data.frame(yall)

### Sum up the counts of methylated and unmethylated reads to get the total read coverage at each CpG site for each sample ###
Methylation <- gl(2, 1, ncol(yall), labels = c("Me", "Un"))
Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]
Coverage <- Me + Un
head(Coverage)

#### No filtering for matrix imputation
n=0
keep_foot <- rowSums(Coverage >= n) >= 0
table(keep_foot)

### DGEList object is subsetted to retain only the filtered loci ###
y_foot <- yall[keep_foot,, keep.lib.sizes=FALSE]

### Normalization - set the library sizes for each sample to be the average of the total read counts for the methylated and unmethylated libraries ###
TotalLibSize <- y_foot$samples$lib.size[Methylation=="Me"] +
  +                 y_foot$samples$lib.size[Methylation=="Un"]
y_foot$samples$lib.size <- rep(TotalLibSize, each=2)
y_foot$samples

### Compute the corresponding methylation summary from the methylated and unmethylated counts ###
Me_foot <- y_foot$counts[, Methylation=="Me"]
Un_foot <- y_foot$counts[, Methylation=="Un"]

### Calculating a methylation proportion matrix
prop_meth_matrix_foot <- Me_foot/(Me_foot+Un_foot)

# Convert to a data frame and make column names syntactically valid
prop_df <- as.data.frame(prop_meth_matrix_foot)
colnames(prop_df) <- make.names(colnames(prop_df))

# Also update the sample names in metadata to match
rownames(meta_data_foot) <- make.names(rownames(meta_data_foot))

# Impute missing values using predictive mean matching
imputed_data <- mice(prop_df, m = 1, method = 'pmm', seed = 123)

# Extract the completed dataset
complete_data <- complete(imputed_data, 1)

###################################################


################ Making Figure 2 ##########
setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/outlier detection/")
pca_scores<-read.csv("all_imp_scores_pca.csv")
meta_data<-read.csv("samples_65.csv")

###GETTING LENGTH DATA FOR ALL 65 SAMPLES

sample_info<-read.csv("full_sample.csv")
all_length<-meta_data %>% 
  left_join(sample_info,
            by=c("meth"="sample.ID",
                 "tissue"="Tissue")) %>% 
  dplyr::select(sample,meth,tissue,origin_site,transplant_site,GROUP_CODE.x,
         LENGTH_INIT..mm.,LENGTH_FINAL..mm.)
# write_csv(all_length,"sample_65_wlength.csv")

##getting the PC contributions
# complete_data <-read.csv("all_complete_data.csv")
# imp_pca_result <- prcomp(t(complete_data))
# summary(imp_pca_result)
#pc1, 0.03749  pc2 0.03026

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
  geom_point(size = 8, alpha = 0.85, stroke = 0.2, color = "black") +
  geom_line(aes(group = meth), color = "black", linetype = 1) +
  scale_shape_manual(
    values = c(
      "lowprotected_protected" = 21,
      "lowprotected_exposed" = 22,
      "lowexposed_exposed" = 24,
      "lowexposed_protected" = 23
    ),
    labels = c(
      "lowprotected_protected" = "Protected_Protected",
      "lowprotected_exposed" = "Protected_Exposed",
      "lowexposed_exposed" = "Exposed_Exposed",
      "lowexposed_protected" = "Exposed_Protected"
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
    x = "PC1 (3.75%)",
    y = "PC2 (3.03%)",
    fill = "Tissue Type",
    shape = "Treatment",
    tag = "A)"  
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
pca_imputated_merged




##### calculating Euclidean distance between tissues ###

## formatting
edistance<-pca_scores %>% 
  select(X, PC1,PC2, Tissue)

edistance$SampleID <- rownames(edistance)
rownames(edistance) <- NULL

edistance$SampleID  <- gsub("^X", "", edistance$SampleID  )

##calculating distance between the same individual

same_indiv_distances <- edistance %>%
  group_by(X) %>%
  filter(dplyr::n() == 2) %>%  
  arrange(X, Tissue) %>%
  summarise(distance = sqrt((PC1[1] - PC1[2])^2 + (PC2[1] - PC2[2])^2)) %>%
  pull(distance)

mean_same_indiv <- mean(same_indiv_distances)
print(mean_same_indiv)


###calculating distance between different individuals

coords <- edistance %>%
  select(X, PC1, PC2)

# making all pairwise combinations
pair_indices <- combn(nrow(coords), 2)

# Compute distances between different tissues of different individuals
between_distances <- apply(pair_indices, 2, function(idx) {
  i <- idx[1]
  j <- idx[2]
  
  if (coords$X[i] != coords$X[j]) {   # only between individuals
    dist <- sqrt((coords$PC1[i] - coords$PC1[j])^2 +
                   (coords$PC2[i] - coords$PC2[j])^2)
    
    # record tissue info for later summaries
    data.frame(
      id1 = coords$X[i],
      id2 = coords$X[j],
      tissue1 = coords$Tissue[i],
      tissue2 = coords$Tissue[j],
      distance = dist
    )
  } else {
    NULL
  }
})

between_df <- do.call(rbind, between_distances)
between_distances <- na.omit(between_df )
mean_diff_indiv <- mean(between_distances$distance)


###### Mantel test ######

# PCA distance matrix across all samples
pca_dist <- dist(edistance[, c("PC1","PC2")])

# 0 if same individual, 1 if different
id_dist <- as.dist( outer(edistance$X, edistance$X,
                          function(a,b) as.integer(a != b)) )

mantel_test <- mantel(pca_dist, id_dist,
                    method = "pearson", permutations = 999)

mantel_test
############

########## regressing PC1 and length

length_data<-read_csv("sample_65_wlength.csv")

#using imp_scores, has all the loading data for 69 or 56 samples
# pca_scores$X<-pca_scores$X%>% 
#   gsub("^X", "", .) %>%        # Remove the initial "X"
#   gsub("_.*$", "", .) %>% 
#   gsub("\\.", "-", .)  

colnames(length_data)[1] <- "sample_name"
PC <- length_data %>% 
  left_join(pca_data, by = c("sample_name"="sample"))

pc1_length<-ggplot(PC, aes(x = `LENGTH_FINAL..mm.`, y = PC1)) + #plotting shell length and PC1 loadings
  geom_point(aes(color= tissue.x), alpha = 0.85, stroke = 0.2,
             size=8) +  # Keep points colored by tissue
  geom_line(aes(group = `meth.x`), color = "black", linetype = 1) + 
  geom_smooth(method = "lm", size = 0.8, alpha = 0.3, color = "blue") +  # Single regression line for all shell length and PC1
  labs(
    y = "PC1", #label y axis
    x = "Final shell length (mm)", #label x axis
    color = "Tissue Type", #specify color of points
    tag = "B)"  
  ) + 
  scale_color_manual(values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e"),
                     labels = c("F" = "Foot", "G" = "Gill"))+ #specify color of points
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    axis.text.x = element_text(size = 20,color="black"), # X-axis text size
    axis.text.y = element_text(size = 20,color="black"), # Y-axis text size
    axis.title.x = element_text(size = 24,color="black"),  # Axis title size
    axis.title.y = element_text(size = 24,color="black"),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14,color="black"), # Legend text size
    legend.title = element_text(size = 16,color="black"),# Legend title size
    panel.grid = element_blank(),         # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none" ,
    legend.box.margin = margin(0, -10, 0, 0)) 

pc1_length


# PC$`sample ID`<- sub("-.*", "",PC$`sample ID`)
# PC_unique <- PC[!duplicated(PC$`sample ID`), ] 

# lm_length_old<-lm(PC1~`LENGTH_FINAL (mm)` + Tissue.x, data=PC)
lm_length<-lmer(PC1~`LENGTH_FINAL..mm.` + tissue.x + (1|`meth.x`), data=PC)
summary(lm_length)

lm_length<-lmer(PC2~`LENGTH_FINAL..mm.` + tissue.x + (1|`meth.x`), data=PC)
summary(lm_length)

# Save both panels for figure 2
p1_fix <- pca_imputated_merged +
  labs(tag = "A.") +
  theme(plot.margin = margin(20,20,20,20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.01, 0.99))

p2_fix <- pc1_length +
  labs(tag = "B.") +
  theme(plot.margin = margin(20,20,20,20),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.01, 0.99))

combined <- (p1_fix | guide_area() | p2_fix) +
  plot_layout(widths = c(1, 0.42, 1), guides = "collect")

# ggsave(
#   "/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/submission/fig2.png",
#   plot = combined, width = 20, height = 8, dpi = 300, bg = "white"
# )

##### Making figure S2 ###
 
pc2_length<-ggplot(PC, aes(x = `LENGTH_FINAL..mm.`, y = PC2)) + #plotting shell length and PC1 loadings
   geom_point(aes(color= tissue.x), alpha = 0.85, stroke = 0.2,
              size=8) +  # Keep points colored by tissue
   geom_line(aes(group = `meth.x`), color = "black", linetype = 1) + 
   geom_smooth(method = "lm", size = 0.8, alpha = 0.3, color = "blue") +  # Single regression line for all shell length and PC1
   labs(
     y = "PC2", #label y axis
     x = "Final shell length (mm)", #label x axis
     color = "Tissue Type" #specify color of points
   ) + 
   scale_color_manual(values = c("Foot" = "#3bbdcc", "Gill" = "#e68d4e"),
                      labels = c("F" = "Foot", "G" = "Gill"))+ #specify color of points
   theme_bw() +
   guides(color = guide_legend(override.aes = list(size = 3))) +
   theme(
     axis.text.x = element_text(size = 20,color="black"), # X-axis text size
     axis.text.y = element_text(size = 20,color="black"), # Y-axis text size
     axis.title.x = element_text(size = 24,color="black"),  # Axis title size
     axis.title.y = element_text(size = 24,color="black"),
     strip.text = element_text(size = 16),
     legend.text = element_text(size = 14,color="black"), # Legend text size
     legend.title = element_text(size = 16,color="black"),# Legend title size
     panel.grid = element_blank(),         # Remove grid lines
     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
     legend.position = "none" ,
     legend.box.margin = margin(0, -10, 0, 0)) 
 
 pc2_length
 

# ggsave("figS2.png", plot = pc2_length, width = 10, height = 8, dpi = 300)
