### Aim of this script: generating average log fold change values for each gene that overlaps with CpGs for use in GOMWU ###

### Foot subset ###

### Load results from DM analyses ###
all_cpg_foot<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/R_code/all_CpG_dm_foot_output.csv")


### reformatting ###
colnames(all_cpg_foot)[colnames(all_cpg_foot) == "...1"] <- "combined"

all_cpg_foot$Chr <- sapply(strsplit(as.character(all_cpg_foot$sample), "-"), `[`, 1)  # Part before '-'
all_cpg_foot$Locus <- sapply(strsplit(as.character(all_cpg_foot$sample), "-"), `[`, 2) 

all_cpg_foot$Locus<-as.numeric(all_cpg_foot$Locus)+1

all_cpg_foot_origin<-all_cpg_foot %>% 
  filter(Treat =="Origin")

all_cpg_foot_transplant<-all_cpg_foot %>% 
  filter(Treat =="Transplant")

### Load gene name information associated with each CpG ###
foot_origin_gene_name<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/GO Enrichment/gene_name_DMsites_foot_origin_corrected.csv")
foot_transplant_gene_name<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/GO Enrichment/gene_name_DMsites_foot_transplant_corrected.csv")

## Getting average LFC for origin effect ###

lfg_gene_foot_origin<-all_cpg_foot_origin %>% 
  left_join(foot_origin_gene_name,
            by=c("Chr","Locus")) %>% 
  dplyr::select(-combined,-sample)

average_lfg_foot_origin<-lfg_gene_foot_origin %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

#write_csv(average_lfg_foot_origin,"average_lfg_foot_origin.csv")
  

## Getting average LFC for transplant effect ###

lfg_gene_foot_transplant<-all_cpg_foot_transplant %>% 
  left_join(foot_transplant_gene_name,
            by=c("Chr","Locus")) %>% 
  dplyr::select(-combined,-sample)

average_lfg_foot_transplant<-lfg_gene_foot_transplant %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

#write_csv(average_lfg_foot_transplant,"average_lfg_foot_transplant.csv")


### Gill subset ###


### Load results from DM analyses ###
all_cpg_gill<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/R_code/all_CpG_dm_gill_output.csv")
colnames(all_cpg_gill)[colnames(all_cpg_gill) == "...1"] <- "combined"

### reformatting ###
all_cpg_gill$Chr <- sapply(strsplit(as.character(all_cpg_gill$sample), "-"), `[`, 1)  # Part before '-'
all_cpg_gill$Locus <- sapply(strsplit(as.character(all_cpg_gill$sample), "-"), `[`, 2) 
all_cpg_gill$Locus<-as.numeric(all_cpg_gill$Locus)+1


all_cpg_gill_origin<-all_cpg_gill %>% 
  filter(Treat =="Origin")

all_cpg_gill_transplant<-all_cpg_gill %>% 
  filter(Treat =="Transplant")

### Load gene name information associated with each CpG ###
gill_origin_gene_name<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/GO Enrichment/gene_name_DMsites_gill_origin_corrected.csv")
gill_transplant_gene_name<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/GO Enrichment/gene_name_DMsites_gill_transplant_corrected.csv")

## Getting average LFC for origin effect ###

lfg_gene_gill_origin<-all_cpg_gill_origin %>% 
  left_join(gill_origin_gene_name,
            by=c("Chr","Locus")) %>% 
  select(-combined,-sample)

average_lfg_gill_origin<-lfg_gene_gill_origin %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

#write_csv(average_lfg_gill_origin,"average_lfg_gill_origin.csv")

## Getting average LFC for transplant effect ###

lfg_gene_gill_transplant<-all_cpg_gill_transplant %>% 
  left_join(gill_transplant_gene_name,
            by=c("Chr","Locus")) %>% 
  select(-combined,-sample)

average_lfg_gill_transplant<-lfg_gene_gill_transplant %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

#write_csv(average_lfg_gill_transplant,"average_lfg_gill_transplant.csv")



### Apply filtering same as those in GOseq analyses ###

go_gill<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/category_list_parent_filtered_gill.csv")
go_foot<-read_csv("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/category_list_parent_filtered_foot.csv")

go_gill_df <- go_gill %>%
  group_by(gene_id) %>%
  summarise(GO = paste(unique(Parent_Term), collapse = ";"), .groups = "drop")


go_foot_df <- go_foot %>%
  group_by(gene_id) %>%
  summarise(GO = paste(unique(Parent_Term), collapse = ";"), .groups = "drop")

### Average LFC after applying filtering threshold ###

lfg_gene_foot_origin_filtered<-lfg_gene_foot_origin %>% 
  filter(gene_id %in% go_foot_df$gene_id)

average_lfg_foot_origin_filtered<-lfg_gene_foot_origin_filtered%>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

lfg_gene_foot_transplant_filtered<-lfg_gene_foot_transplant %>% 
  filter(gene_id %in% go_foot_df$gene_id)

average_lfg_foot_transplant_filtered<-lfg_gene_foot_transplant_filtered%>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

#write_csv(average_lfg_foot_origin_filtered,"/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/average_lfg_foot_origin_filtered.csv")
#write_csv(average_lfg_foot_transplant_filtered,"/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/average_lfg_foot_transplant_filtered.csv")
# 

lfg_gene_gill_transplant_filtered <- lfg_gene_gill_transplant %>%
  filter(gene_id %in% go_gill_df$gene_id)

average_lfg_gill_transplant_filtered<-lfg_gene_gill_transplant_filtered %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

lfg_gene_gill_origin_filtered <-lfg_gene_gill_origin%>% 
  filter(gene_id %in% go_gill_df$gene_id)

average_lfg_gill_origin_filtered <- lfg_gene_gill_origin_filtered %>% 
  group_by(gene_id) %>% 
  summarize(mean=mean(logFC)) %>% 
  na.omit()

# write_csv(average_lfg_gill_transplant_filtered,"/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/average_lfg_gill_transplant_filtered.csv")
# write_csv(average_lfg_gill_origin_filtered,"/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/MethylMussel_GO/GO_MWU/average_lfg_gill_origin_filtered.csv")

