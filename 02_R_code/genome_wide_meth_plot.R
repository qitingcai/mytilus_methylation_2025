###Aim of this script: test whether genome wide DNA methylation proportion differ by tissue type or treatment group, plot the dara ###


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
