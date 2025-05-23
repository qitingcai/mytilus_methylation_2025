### Aim of this script: make bar plots that shows the proportion of CpGs or DM CpGs that fall into each genomic feature ###

### plot showing proportions of DM CpGs falling into each genomic features ###
dm_bar_plot<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/writing/dm_genic_feature.csv")

dm_bar_plot<- dm_bar_plot %>% 
  mutate(nonDM=Total-DM_count)

# Define custom colors (feature = color)
custom_colors <- c(
  "Intron" = "#0662AE",
  "Exon" = "#ABCDE9",
  "Promoter" = "#1F3161",
  "Intergenic" = "#3792D0"
)

library("ggthemes")

library(dplyr)

# Create a new row for the "Total" feature in each group
total_bars <- dm_bar_plot %>%
  group_by(Type, tissue) %>%
  summarise(
    DM_count = sum(DM_count),
    nonDM_count = sum(nonDM_count),
    Total = sum(Total),
    .groups = "drop"
  ) %>%
  mutate(
    feature = "Total",
    DM_proportion = DM_count / Total,
    non_DM_proportion = nonDM_count / Total
  )



# # Plot the data
ggplot(dm_bar_plot, aes(x = Type, y = DM_proportion, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~tissue) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +  # Ensures all colors show
  labs(
    title = "",
    y = "Proportion (%)",
    x = "",
    fill = "Feature"
  ) +
  theme_few()




foot<-dm_bar_plot %>% 
  filter(tissue=="Foot") %>% 
  pivot_longer(
    cols = c(DM_proportion, Total),
    names_to = "DM_status",
    values_to = "proportion"
  ) %>%
  mutate(
    DM_status = recode(DM_status,
                       "DM_proportion" = "DM",
                       "Total" = "Total")
  )


p1 <- ggplot(foot, aes(x = DM_status, y = proportion, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Type) +
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
    axis.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size =25),
    strip.text = element_text(color = "black", size = 25),
    plot.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p1


gill<-dm_bar_plot %>% 
  filter(tissue=="Gill") %>% 
  pivot_longer(
    cols = c(DM_proportion, Total),
    names_to = "DM_status",
    values_to = "proportion"
  ) %>%
  mutate(
    DM_status = recode(DM_status,
                       "DM_proportion" = "DM",
                       "Total" = "Total")
  )

p2<-ggplot(gill, aes(x = DM_status, y = proportion, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Type) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +  # Ensures all colors show
  labs(
    title = "",
    y = "",
    x = "",
    fill = "Feature"
  ) +
  theme_few() +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right",
    axis.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 25),
    strip.text = element_text(color = "black", size = 25),
    plot.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
p2

# p1 <- p1 +
#   theme(
#     plot.margin = margin(10, 10, 10, 10),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14)
#   )
# 
# p2 <- p2 +
#   theme(
#     plot.margin = margin(10, 10, 10, 10),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14)
#   )


# library(patchwork)
# 
# foot_DM <- p1 +plot_annotation(title = "A. Foot",
#                                     theme = theme(plot.title = element_text(face= "bold", size = 20)))
# 
# gill_DM <- p2 +
#   plot_annotation(title = "B. Gill",
#                   theme = theme(plot.title = element_text(face= "bold", size = 20)))
# 
# 
# 
# combined <- wrap_elements(foot_DM) + wrap_elements(gill_DM) +
#   plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
#   plot_annotation(tag_levels = NULL)  # turn off A/B tags
# 
# combined
# 
# # Reapply shared layout and theme
# shared_theme <- theme_few() +
#   theme(
#     axis.text = element_text(color = "black", size = 12),
#     axis.title = element_text(color = "black", size = 15),
#     strip.text = element_text(color = "black", size = 15),
#     plot.title = element_text(face = "bold", size = 20, hjust = 0),  # left-align title
#     axis.ticks = element_line(color = "black"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     legend.text = element_text(size = 13),
#     legend.title = element_text(size = 15),
#     legend.position = "right"
#   )
# 
# # Apply titles as ggtitle
# p1 <- p1 +
#   ggtitle("A. Foot") +
#   shared_theme +
#   labs(y = "Proportion", x = "", fill = "Feature")
# 
# p2 <- p2 +
#   ggtitle("B. Gill") +
#   shared_theme +
#   labs(y = "", x = "", fill = "Feature")
# 
# # Combine with shared legend and equal width
# combined <- p1 + p2 +
#   plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
#   plot_annotation(tag_levels = NULL)
# 
# combined



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
    plot.margin = margin(10, 10, 10, 2)
  )

# Combine with cowplot
library(cowplot)

plots_combined <- plot_grid(p1_clean, p2_clean, ncol = 2, rel_widths = c(1, 1), align = "hv")

# Combine plots with smaller relative widths
plots_combined <- plot_grid(
  p1_clean, p2_clean,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "hv"
)

combined<-ggdraw() +
  draw_plot(plots_combined, 0, 0, 1, 0.95) +
  draw_plot_label(
    label = c("A. Foot", "B. Gill"),
    x = c(0.025, 0.528),  # ← this is the key tweak
    y = c(1, 1),
    hjust = 0,
    fontface = "bold",
    size = 26
  )
combined

ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/submission/DM_CpG_combined.png", combined, device = "png", width = 20, height = 11, dpi = 300, units = "in")
ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/submission/DM_CpG_combined.pdf", combined, device = "pdf", width = 14, height = 10, dpi = 300, units = "in")





#chisquare
#to look at overall distribution and whether the patterns differ between tissue
table_data <- xtabs(DM_count ~ feature + tissue, data = dm_bar_plot)
fisher.test(table_data)

#test whether gene body and promoter greater in foot vs. gill

gene_body_promoter <- colSums(table_data[c("exon", "intron", "promoter"), ])
intergenic <- table_data["intergenic", ]

# Create 2x2 matrix
collapsed_table <- rbind(GeneBody_Promoter = gene_body_promoter,
                         Intergenic = intergenic)

fisher.test(collapsed_table)



#look at dm cpgs compared to all available cpgs
gill_trans <- dm_bar_plot %>% filter(tissue == "Gill", Type == "Transplant")
gill_origin <- dm_bar_plot %>% filter(tissue == "Gill", Type == "Origin")


foot_trans <- dm_bar_plot %>% filter(tissue == "Foot", Type == "Transplant")
foot_origin <- dm_bar_plot %>% filter(tissue == "Foot", Type == "Origin")

# For each feature, compare to rest
enrichment_results <- foot_origin%>%
  rowwise() %>%
  mutate(
    other_DM = sum(foot_origin$DM_count) - DM_count, #dm found in other features
    other_total = sum(foot_origin$Total) - Total, #other available features
     p_value = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$p.value,
     odds_ratio = fisher.test(matrix(c(DM_count, Total, other_DM, other_total), nrow = 2))$estimate
  ) %>%
  ungroup() %>% 
  select(feature, DM_count, Total, p_value, odds_ratio)


enrichment_results <- enrichment_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

enrichment_results

# # ###GLM with negative binomial distribution to see if count in each category differ
# # library(MASS)
#  proportion_glm<-MASS::glm.nb(count~feature*tissue, data=dm_bar_plot)
#  summary(proportion_glm)
#  library(emmeans)
# emmeans(proportion_glm, pairwise ~ feature*tissue)
#  emmeans(proportion_glm, pairwise ~ tissue)
# # contrast(emm, method = "pairwise", by = "feature", adjust = "none")
# # 
# # emmeans(proportion_glm, pairwise ~ type*tissue)



#include 3 bars instead of 4

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
                             "Origin_DM_count" = "Origin DM",
                             "Transplant_DM_count" = "Transplant DM",
                             "Total" = "Total CpGs")
  )


dm_long <- dm_long %>%
  mutate(
    Type_Categories = factor(Type_Categories, levels = c("Origin DM", "Transplant DM", "Total CpGs"))
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
    axis.text = element_text(color = "black", size = 23),
    axis.title = element_text(color = "black", size =23),
    strip.text = element_text(color = "black", size = 25),
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
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22.8),
    axis.text = element_text(color = "black", size = 22),
    axis.title = element_text(color = "black", size =23),
    strip.text = element_text(color = "black", size = 25),
    plot.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
p2




# # Plot the data
ggplot(dm_long, aes(x = Type_Categories, y = Count, fill = feature)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~tissue) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +  # Ensures all colors show
  labs(
    title = "",
    y = "Proportion (%)",
    x = "",
    fill = "Feature"
  ) +
  theme_few()+ theme(
    strip.text = element_text(face = "bold", size = 12)
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
    plot.margin = margin(10, 10, 10, 2)
  )

# Combine with cowplot
library(cowplot)

plots_combined <- plot_grid(p1_clean, p2_clean, ncol = 2, rel_widths = c(1, 1), align = "hv")

# Combine plots with smaller relative widths
plots_combined <- plot_grid(
  p1_clean, p2_clean,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "hv"
)

combined<-ggdraw() +
  draw_plot(plots_combined, 0, 0, 1, 0.95) +
  draw_plot_label(
    label = c("A. Foot", "B. Gill"),
    x = c(0.008, 0.48),  # ← this is the key tweak
    y = c(1, 1),
    hjust = 0,
    fontface = "bold",
    size = 26
  )
combined


ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/submission/DM_CpG_combined_fix0513.png", combined, device = "png", width = 20, height = 11, dpi = 300, units = "in")
ggsave("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/submission/DM_CpG_combined_fix0513.pdf", combined, device = "pdf", width = 14, height = 10, dpi = 300, units = "in")

