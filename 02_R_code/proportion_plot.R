### Aim of this script: make bar plots that shows the proportion of CpGs or DM CpGs that fall into each genomic feature ###

### plot showing proportions of DM CpGs falling into each genomic features ###
dm_bar_plot<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/manuscript/writing/dm_genic_feature.csv")

# Define custom colors (feature = color)
custom_colors <- c(
  "intron" = "#0662AE",
  "exon" = "#ABCDE9",
  "promoter" = "#1F3161",
  "intergenic" = "#3792D0"
)

library("ggthemes")

# Plot the data
ggplot(dm_bar_plot, aes(x = type, y = proportion, fill = feature)) +
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
