### Aim of this script: make bar plots that shows the proportion of CpGs or DM CpGs that fall into each genomic feature ###

### plot showing proportions of DM CpGs falling into each genomic features ###
##updated file with interaction term rerun

setwd("/Users/qcai/mytilus_github/Data")

dm_bar_plot<-read_csv("fisher_test_input_final_revised.csv")

dm_bar_plot<- dm_bar_plot %>% 
  mutate(nonDM=Total-DM_count)

# Define custom colors (feature = color)
custom_colors <- c(
  "Intron" = "#0662AE",
  "Exon" = "#ABCDE9",
  "Promoter" = "#1F3161",
  "Intergenic" = "#3792D0",
  "UTR"="#d8e0e6"
)

library("ggthemes")
library(dplyr)
library(dplyr)

#look at dm cpgs compared to background
gill_trans <- dm_bar_plot %>% filter(tissue == "Gill", Type == "Transplant")
gill_origin <- dm_bar_plot %>% filter(tissue == "Gill", Type == "Origin")
# gill_interact <- dm_bar_plot %>% filter(tissue == "Gill", Type == "Interacton")

foot_trans <- dm_bar_plot %>% filter(tissue == "Foot", Type == "Transplant")
foot_origin <- dm_bar_plot %>% filter(tissue == "Foot", Type == "Origin")
# foot_interact <- dm_bar_plot %>% filter(tissue == "Foot", Type == "Interacton")

origin_tissues<- dm_bar_plot %>% 
  filter(Type =="Origin")
transplant_tissues<- dm_bar_plot %>% 
  filter(Type =="Transplant")


#to look at overall distribution and whether the patterns differ between tissue combined
table_data_origin <- xtabs(DM_count ~ feature + tissue, data = origin_tissues)
fisher.test(table_data_origin)

table_data_trans <- xtabs(DM_count ~ feature + tissue, data = transplant_tissues)
fisher.test(table_data_trans)


# Transplant site comparing distribution of each feature vs. all other features in foot and gills

features_trans <- rownames(table_data_trans)
results <- lapply(features_trans, function(f) {
  # rows to use for "all other features"
  others <- setdiff(features_trans, f)
  
  feature_vs_rest <- rbind(
    Feature    = table_data_trans[f, , drop = FALSE],
    All_others = colSums(table_data_trans[others, , drop = FALSE])
  )
  
  test <- fisher.test(feature_vs_rest)  # works if there are exactly 2 columns
  
  data.frame(
    feature    = f,
    odds_ratio = if (!is.null(test$estimate)) unname(test$estimate) else NA_real_,
    p_value    = test$p.value,
    row.names  = NULL
  )
})
results_df <- bind_rows(results)
results_df <- results_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

results_df




# Origin site comparing distribution of each feature vs. all other features in foot and gills

features_origin <- rownames(table_data_origin)
results <- lapply(features_origin , function(f) {
  # rows to use for "all other features"
  others <- setdiff(features_origin, f)
  
  feature_vs_rest <- rbind(
    Feature    = table_data_origin[f, , drop = FALSE],
    All_others = colSums(table_data_origin[others, , drop = FALSE])
  )
  
  test <- fisher.test(feature_vs_rest)  # works if there are exactly 2 columns
  
  data.frame(
    feature    = f,
    odds_ratio = if (!is.null(test$estimate)) unname(test$estimate) else NA_real_,
    p_value    = test$p.value,
    row.names  = NULL
  )
})

results_df <- bind_rows(results)
results_df <- results_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
results_df



# For each feature, compare to rest
enrichment_results <- foot_origin%>%
  rowwise() %>%
  mutate(
    other_DM = sum(foot_origin$DM_count) - DM_count, #dm found in other features
    other_total = sum(foot_origin$Total) - Total, #other available features
    p_value = fisher.test(matrix(c(DM_count, Total , other_DM, other_total), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$estimate
  ) %>%
  ungroup() %>% 
  dplyr::select(feature, DM_count, Total, p_value, odds_ratio)


enrichment_results <- enrichment_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

enrichment_results



# For each feature, compare to rest
enrichment_results <- foot_trans%>%
  rowwise() %>%
  mutate(
    other_DM = sum(foot_trans$DM_count) - DM_count, #dm found in other features
    other_total = sum(foot_trans$Total) - Total, #other available features
    p_value = fisher.test(matrix(c(DM_count, Total , other_DM, other_total), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$estimate
  ) %>%
  ungroup() %>% 
  dplyr::select(feature, DM_count, Total, p_value, odds_ratio)


enrichment_results <- enrichment_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

enrichment_results




# For each feature, compare to rest
enrichment_results <- gill_origin%>%
  rowwise() %>%
  mutate(
    other_DM = sum(gill_origin$DM_count) - DM_count, #dm found in other features
    other_total = sum(gill_origin$Total) - Total, #other available features
    p_value = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$estimate
  ) %>%
  ungroup() %>% 
  dplyr::select(feature, DM_count, Total, p_value, odds_ratio)


enrichment_results <- enrichment_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

enrichment_results



# For each feature, compare to rest
enrichment_results <- gill_trans%>%
  rowwise() %>%
  mutate(
    other_DM = sum(gill_trans$DM_count) - DM_count, #dm found in other features
    other_total = sum(gill_trans$Total) - Total, #other available features
    p_value = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$estimate
  ) %>%
  ungroup() %>% 
  dplyr::select(feature, DM_count, Total, p_value, odds_ratio)


enrichment_results <- enrichment_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

enrichment_results



# Pivot DM_count, nonDM_count, and Total into long format
library(dplyr)
library(tidyr)


dm_long <- dm_bar_plot %>%
  dplyr::select(feature, tissue, Type, DM_count, nonDM_count, Total) %>%
  pivot_longer(
    cols = c(DM_count, nonDM_count, Total),
    names_to = "Categories",
    values_to = "Count"
  ) %>%
  # Create a new Type_Categories column, only combining for DM_count
  mutate(
    Type_Categories = ifelse(
      Categories == "Total",
      "Total",  # keep as "Total"
      paste(Type, Categories, sep = "_")  # e.g., "Origin_DM_count"
    )
  ) %>%
  # Filter out nonDM if desired
  filter(Categories != "nonDM_count")


dm_long <- dm_long %>%
  mutate(
    Type_Categories = recode(Type_Categories,
                             "Origin_DM_count" = "Origin",
                             "Transplant_DM_count" = "Transplant",
                             "Interacton_DM_count" = "Interaction",
                             "Total" = "Total"))


dm_long <- dm_long %>%
  mutate(
    Type_Categories = factor(Type_Categories, levels = c("Origin", "Transplant", "Interaction", "Total"))
  )


foot<-dm_long %>% 
  filter(tissue=="Foot") 

gill<-dm_long %>% 
  filter(tissue=="Gill") 


p1 <- ggplot(foot, aes(x = Type_Categories, y = Count, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  labs(
    title = "",
    y = "Proportion",
    x = "",
    fill = "Feature"
  ) +
  theme_few() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(color = "black", size = 12),  
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size =20),
    strip.text = element_text(color = "black", size = 20),
    plot.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p1


p2 <- ggplot(gill, aes(x = Type_Categories, y = Count, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  labs(
    title = "",
    y = "Proportion",
    x = "",
    fill = "Feature"
  ) +
  theme_few() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(color = "black", size = 12),  
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size =20),
    strip.text = element_text(color = "black", size = 20),
    plot.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p2

p2<-p2 + annotate(
  "text",
  x = "Origin",
  y = 0.89,
  label = "*",
  size = 10,
  color = "black"
) + annotate(
  "text",
  x = "Origin",
  y = 0.62,
  label = "*",
  size = 10,
  color = "black"
) + annotate(
  "text",
  x = "Transplant",
  y = 0.93,
  label = "*",
  size = 10,
  color = "black"
)  + annotate(
  "text",
  x = "Transplant",
  y = 0.535,
  label = "*",
  size = 10,
  color = "black"
)


# Add margin to ensure y-axis label is visible
p1_clean <- p1 +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(10, 1, 10, 0.2),
    legend.position = "none"# Top, Right, Bottom, Left
  )

p2_clean <- p2 +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(10, 2, 10, 0.2)
  )


# Combine with cowplot
library(cowplot)

plots_combined <- plot_grid(p1_clean, p2_clean, ncol = 2, rel_widths = c(1, 1), align = "hv")

# # Combine plots with smaller relative widths
plots_combined <- plot_grid(
  p1_clean, p2_clean,
  ncol = 2,
  align = "hv",
  axis = "tb" 
)

combined<-ggdraw(ggdraw(clip = "off") ) +
  theme(plot.background = element_rect(fill = "white", colour = NA))+
  draw_plot(plots_combined, 0, 0, 1, 0.95) +
  draw_plot_label(
    label = c("A. Foot", "B. Gill"),
    x = c(0.008, 0.48),
    y = c(0.999, 0.999),
    hjust = 0,
    fontface = "bold",
    size = 20
  )



combined

# ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/DM_CpG_combined.png", combined, device = "png", width = 12, height = 6, dpi = 300, units = "in")
ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/Round1_revision/DM_CpG_combined_revised.png", combined, device = "png", width = 12, height = 6, dpi = 300, units = "in")


#### glm between tissues ####

#### see if log fold changes differ between treatment and tissues

library(emmeans)
merged_all_DM_gill<-read_csv("/Users/qcai/mytilus_github/Data/merged_all_DM_gill.csv") %>% mutate(tissue ="Gill")
merged_all_DM_foot<-read_csv("/Users/qcai/mytilus_github/Data/merged_all_DM_foot.csv") %>% mutate(tissue ="foot")


merged_unique_foot <- merged_all_DM_foot %>%
  distinct(combined,Treat, .keep_all = TRUE)

merged_unique_gill <- merged_all_DM_gill %>%
  distinct(combined,Treat, .keep_all = TRUE)


lfc_tissues<-rbind(merged_unique_gill,merged_unique_foot)

fit <- lm(transplant_logFC ~ tissue * Treat, data = lfc_tissues)
summary(fit)

emmeans(fit, pairwise ~ Treat | tissue)   # treatment means within each tissue
emmeans(fit, pairwise ~  tissue | Treat) 


abs_fit <- lm(abs(transplant_logFC) ~ tissue * Treat , data = lfc_tissues)
summary(abs_fit)
emmeans(abs_fit, pairwise ~ Treat | tissue)   # treatment means within each tissue
emmeans(abs_fit, pairwise ~  tissue | Treat) 

