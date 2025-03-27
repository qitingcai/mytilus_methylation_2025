### foot specific ###
############ Differential Methylation Analysis ############

### Load libraries ###
library(vegan)
library(edgeR)
library(tidyverse)
library(SYNCSA)
library(mice)
library(ape)
library(rtracklayer)
library(genomation)
library(plyranges)
library(GenomicRanges)
library(lme4)

### Set working directory, include all (top5 for each treatments) foot coverage files ###
setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/tissue_specific/foot")

### Load metadata file that include sample names, load all coverage files into R ###
meta_data_foot<-read.delim("foot_metadata.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_foot <- row.names(meta_data_foot)
files_foot <- paste0(Sample_foot,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

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

### Filtering to only include samples with a coverage of at least 3, for 66% of samples, total of 14 samples ###
n=3
keep_foot <- rowSums(Coverage >= n) >= 14
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

### We want to use the exposed site as reference level for both transplant and origin site effect analyses, need to relevel for transplant site effect ###
meta_data_foot$transplant_site<- relevel(factor(meta_data_foot$transplant_site), ref = "exposed")

### Create a design matrix, using origin site, transplant site, and final shell length as fixed effects ###
designSL_foot <- model.matrix(~0+ origin_site + transplant_site +
                           LENGTH_FINAL..mm., 
                         data=meta_data_foot)

### Expand to the full design matrix modeling the sample and methylation effects ###
design_foot <- modelMatrixMeth(designSL_foot)

### Dispersion estimation ###
y_foot <- estimateDisp(y_foot, design = design_foot, robust = TRUE)

### Create the BCV plot ### 
plotBCV(y_foot)

### Testing for differentially methylated CpG loci ###

### fit NB GLMs for all the CpG loci using the glmFit function in edgeR ###
fit_foot <- glmFit(y_foot, design_foot)

### Testing for differentially methylated CpG sites between different treatment groups using likelihood ratio tests ###
### Origin site effects ###
contr_origin_foot <- makeContrasts(Origin = origin_siteprotected-origin_siteexposed, levels = design_foot)
lrt_origin_foot <- glmLRT(fit_foot, contrast=contr_origin_foot)

### Transplant site effects (exposed is the reference level) ###
contr_trans_foot <- makeContrasts(Transplant = transplant_siteprotected, levels = design_foot)
lrt_trans_foot <- glmLRT(fit_foot, contrast=contr_trans_foot)

### Wrangle data for generating volano plot of CpG DM ###
### Correct p-values using BH method ###
lrt_origin_foot$table$FDR <- p.adjust( lrt_origin_foot$table$PValue, method = "BH" )
lrt_trans_foot$table$FDR <- p.adjust( lrt_trans_foot$table$PValue, method = "BH" )

### Apply logical variable for significant DM ###
lrt_origin_foot$table$Sig <- ifelse(lrt_origin_foot$table$FDR< 0.05, TRUE, FALSE )
lrt_trans_foot$table$Sig <- ifelse(lrt_trans_foot$table$FDR< 0.05, TRUE, FALSE )

### Apply binary variable for hyper/hypomethylation or "Up" vs "Down" ###
lrt_origin_foot$table$Dir <- ifelse(lrt_origin_foot$table$logFC > 0, "Up", "Down" )
lrt_trans_foot$table$Dir <- ifelse( lrt_trans_foot$table$logFC > 0, "Up", "Down" )

### Create combined term for significance and fold-change direction of diff meth ###
lrt_origin_foot$table$Sig_Dir <- paste( lrt_origin_foot$table$Sig,
                             lrt_origin_foot$table$Dir,
                             sep = "_" )
lrt_trans_foot$table$Sig_Dir <- paste( lrt_trans_foot$table$Sig,
                             lrt_trans_foot$table$Dir,
                             sep = "_" )

### Create CpG ID and treatment variables ###
lrt_origin_foot$table$Treat <- "Origin"
lrt_trans_foot$table$Treat <- "Transplant"

lrt_origin_foot$table$sample <- rownames(lrt_origin_foot$table)
lrt_trans_foot$table$sample <- rownames(lrt_trans_foot$table)

### Merge origin site and transplant site coefficients ###
all_CpG_dm <- rbind( lrt_origin_foot$table,
                     lrt_trans_foot$table )

### Saving DM sites for transplant site effect and origin site effect as new objects ###
DM_Trans_foot<-all_CpG_dm %>% filter(Sig =="TRUE",Treat =="Transplant")
DM_Origin_foot<-all_CpG_dm %>% filter(Sig =="TRUE", Treat =="Origin")

### Origin site and transplant site associated DM CpG volcano plots ###
ggplot( data = all_CpG_dm, aes(y = -log( as.numeric( FDR ) ), x = as.numeric( logFC ),
                               color = Sig_Dir ) ) +
    geom_point() +
    theme_classic( base_size = 40 ) +
    theme( legend.position = "none",
           strip.background = element_blank() ) +
    scale_color_manual( values = c( "black", "black", "blue", "red" ) ) +
    facet_grid( Treat ~ . , scale = "free") +
    labs( x = "Foot Diff meth", y = "-log FDR" )

#write.csv(all_CpG_dm, "all_CpG_dm_foot_output.csv")

### Find genomics regions that each (DM) CpG falls into ###
### Read in intron annotated gff, see Preprocess_00/add_intron.sh for code for intron annotation ###
gff_intron<-read.gff("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/new_genomic_intron.gff")

### Convert the gff object to a GRanges object using the gffToGRanges() function of the R package genomation ###
all_GRange<-gffToGRanges("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/TOP_5/new_genomic_intron.gff", filter = NULL, zero.based = FALSE, ensembl = FALSE) 


### Divide the annotation GRanges object into different objects filtered for exons and introns ###
exon_GRange<-all_GRange %>% 
  filter(type == c
         ("exon"))
intron_GRange<-all_GRange %>% 
  filter(type == c
         ("intron")) 
         
### Code to find PROMOTER REGIONS: defined as 1kb downstream or upstream of first exons of all genes ###

### Separate positive and negative strand of all exons ###
exon_pos_strand <- exon_GRange %>% filter(strand == "+")
exon_neg_strand <- exon_GRange %>% filter(strand == "-")

### Extracting all first exons of each gene ###
### For positive strands, get the first exons by extracting the smallest start position with slice_min ###
first_exon_pos <- as.data.frame(exon_pos_strand) %>% 
  group_by(gene) %>% 
  slice_min(order_by = start, n = 1, with_ties =FALSE) # Only return one when there are overlapping first exons to avoid duplicates

### For negative strands, extract the first exons using the largest start position with slice_max, (going from right to left) ###
first_exon_neg <- as.data.frame(exon_neg_strand) %>% 
  group_by(gene) %>% 
  slice_max(order_by = start, n = 1, with_ties =FALSE)

### Combine both positive and negative strand first exons ###
first_exon_GRange <- bind_rows(first_exon_pos, first_exon_neg)

### Remove potential duplicated lines ###
first_exon_GRange <- first_exon_GRange %>% distinct()

### Create a dataframe ###
first_exon_GRange_df<-as.data.frame(first_exon_GRange)

### Create new GRange object with only first exons ###
first_exon_GRange <- GRanges(
  seqnames = Rle(first_exon_GRange$seqnames),
  ranges = IRanges(start = first_exon_GRange$start, end = first_exon_GRange$end),
  strand = first_exon_GRange$strand,
  gene_id =first_exon_GRange$gene # genes assigned to promoters are those that are associated with the first exons
)

### For some first exons start site <1000, 1kb upstream returns negative values, adjust length that the exons with start position less than 1000 has a start site of 1  

### Adjust start positions for promoter regions based on strand ###
adjusted_start <- ifelse(strand(first_exon_GRange) == "+",
                         start(first_exon_GRange) - 1001,  # Positive strand: promoter starts from upstream 1kb, not counting the first position of exon
                         end(first_exon_GRange)+1)    # Negative strand: end of first exon (from left to right), add 1 to avoid overlap with exon

### Ensure the adjusted start positions are not less than 1 ###
adjusted_start <- pmax(adjusted_start, 1)

# Adjust end positions for promoters based on strand and start conditions
adjusted_end <- ifelse(
  strand(first_exon_GRange) == "+" & start(first_exon_GRange) == 1,
  1000,  # If positive strand and start is 1, end at 1000
  ifelse(
    strand(first_exon_GRange) == "+",
    start(first_exon_GRange) - 1,  # For positive strand, end 1 bp before the first exon to avoid overlap
    end(first_exon_GRange)+1001  # For negative strand, end 1000 downstream of first exon
  )
)

### Create the GRanges object for promoter regions ###
promoters_Grange<- GRanges(
  seqnames = seqnames(first_exon_GRange),
  ranges = IRanges(
    start = adjusted_start,
    end = adjusted_end
  ),
  strand = strand(first_exon_GRange),
  gene = mcols(first_exon_GRange)$gene_id 
)

### Save promoter GRange as new dataframe ###
promoters_Grange_df<-as.data.frame(promoters_Grange)


### Convert a data frame of coordinates of CpGs to a GRange object using the function makeGRangesFromDataFrame() of the R package GenomicRanges ###

#### Use prop_meth_matrix from above, extract first rows as all the CpG positions across samples ###
prop_meth_matrix_foot <- as.data.frame(prop_meth_matrix_foot)

### Reformating ###
prop_meth_matrix_foot_new <- prop_meth_matrix_foot %>%
  rownames_to_column(var = "name") #change row names to column called name

CpGs_foot <- prop_meth_matrix_foot_new  %>%
  separate(col = 1, into = c("seqid", "start"), sep = "-") #separating chromosome, start site information

CpGs_foot$start<-as.numeric(CpGs_foot$start)

CpGs_foot<-CpGs_foot %>% 
  mutate(start= start+1,
         end =start+1)  #changing from 0 to 1 base format because the annotation file is 1 base, but the coverage files are 0 base format.

### Reorder columns to place 'end' in the third column for easier viewing ###
CpGs_foot <- CpGs_foot[, c(1, 2, ncol(CpGs_foot), 3:(ncol(CpGs_foot)-1))]

### Make cpg a GRange object for all CpG sites ###
CpGs_foot<-makeGRangesFromDataFrame(CpGs_foot)
df_cpg_foot<-as.data.frame(CpGs_foot) ### save a new dataframe for all CpGs ###

### Subsetting CpG GRange objects according to genic features ###
### Exon subset ###
exon_subset_foot<-subsetByOverlaps(CpGs_foot, exon_GRange)
#exon_subset<-as.data.frame(exon_subset) 

### Adding gene id to each entry ###
exon_gff <- subset(all_GRange, type == "exon") #need to be GRange object to subset
exon_id_overlaps_foot <- findOverlaps(exon_subset_foot, exon_gff)

### Extract gene_id for matching exons from the GFF file ###
gene_ids_foot <- mcols(exon_gff)$gene[subjectHits(exon_id_overlaps_foot)]
query_hits_foot <- queryHits(exon_id_overlaps_foot)
# Add gene_id only to those rows in exon_subset that have corresponding gene overlaps

### Initialize with NA for those rows that have no gene overlap ###
exon_subset_gene_id_foot <- rep(NA, length(exon_subset_foot))
# Assign the gene_ids where there is a match
exon_subset_gene_id_foot[query_hits_foot] <- gene_ids_foot

### Add gene_id as a metadata column to exon_subset df ###
mcols(exon_subset_foot)$gene_id <- exon_subset_gene_id_foot
exon_subset_gene_id_df_foot<-as.data.frame(exon_subset_foot)
exon_subset_gene_id_df_foot$combined <- paste(exon_subset_gene_id_df_foot$seqnames, exon_subset_gene_id_df_foot$start, sep = "-") #reformat

### Intron subset ###
intron_subset_foot<-subsetByOverlaps(CpGs_foot, intron_GRange)
#intron_subset_foot<-as.data.frame(intron_subset_foot) 

### Adding gene id ###
intron_gff <- subset(all_GRange, type == "intron") #need to be GRange object to subset
mcols(intron_gff)
intron_id_overlaps_foot <- findOverlaps(intron_subset_foot, intron_gff)

### Extract gene_id for matching exons from the GFF file ###
gene_ids_foot<- mcols(intron_gff)$gene[subjectHits(intron_id_overlaps_foot)]
query_hits <- queryHits(intron_id_overlaps_foot)

### Add gene_id only to those rows in exon_subset that have corresponding overlaps ###
### Initialize with NA for those rows that have no overlap ###
intron_subset_gene_id_foot <- rep(NA, length(intron_subset_foot))

### Assign the gene_ids where there is a match ###
intron_subset_gene_id_foot[query_hits] <- gene_ids_foot

### Add gene_id as a metadata column to intron_subset df ###
mcols(intron_subset_foot)$gene_id <- intron_subset_gene_id_foot
intron_subset_gene_id_df_foot<-as.data.frame(intron_subset_foot)
intron_subset_gene_id_df_foot$combined <- paste(intron_subset_gene_id_df_foot$seqnames, intron_subset_gene_id_df_foot$start, sep = "-") #reformat

### Find promoter overlaps with CpG using the adjusted ranges ###
promoter_GRange_foot <- findOverlaps(CpGs_foot, promoters_Grange)

promoter_subset_foot<-subsetByOverlaps(CpGs_foot, promoters_Grange)
promoter_subset_foot_df<-as.data.frame(promoter_subset_foot) 
#1729 observation

### Adding gene id ###
promoter_gff <- promoters_Grange # need to be GRange object to subset
promoter_id_overlaps_foot <- findOverlaps(promoter_subset_foot, promoter_gff)

### Extract gene_id for matching promoter region (genes associated with first exons, see above, from the GFF file ###
gene_ids_foot<- mcols(promoter_gff)$gene[subjectHits(promoter_id_overlaps_foot)]
query_hits <- queryHits(promoter_id_overlaps_foot)

### Add gene_id only to those rows in promoter_subset_foot that have corresponding overlaps ###
# Initialize with NA for those rows that have no overlap
promoter_subset_gene_id_foot <- rep(NA, length(promoter_subset_foot))

### Assign the gene_ids where there is a match ###
promoter_subset_gene_id_foot[query_hits] <- gene_ids_foot

### Add gene_id as a metadata column to promoter subset df ###
mcols(promoter_subset_foot)$gene_id <- promoter_subset_gene_id_foot
promoter_subset_gene_id_df_foot<-as.data.frame(promoter_subset_foot)
promoter_subset_gene_id_df_foot$combined <- paste(promoter_subset_gene_id_df_foot$seqnames, promoter_subset_gene_id_df_foot$start, sep = "-") #reformat


### FINDING INTERGENIC REGIONS ###
### Filter for all genic regions (exon, intron, UTRs) ###
non_interg_GRange <- subset(all_GRange, type %in% c("exon", "intron", "five_prime_UTR", "three_prime_UTR"))

### Find overlaps between CpG and non-intergenic regions (genic regions) ###
noninterg_overlap <- findOverlaps(CpGs_foot, non_interg_GRange)

### Get the indices of the CpG regions that overlap with genic regions ###
overlapping_indices <- queryHits(noninterg_overlap)

### Remove these overlapping CpG regions to get intergenic regions ###
interg_subset <- CpGs_foot[-overlapping_indices]

### Convert to data frame and make a Grange Object, genes are NAs ###
interg_subset_df_foot <- as.data.frame(interg_subset) %>% mutate(gene_id =NA)
interg_subset_df_foot$combined <- paste(interg_subset_df_foot$seqnames, interg_subset_df_foot$start, sep = "-")
interg_subset_df_foot_Grange<-makeGRangesFromDataFrame(interg_subset_df_foot)


############ Running GLM to assess the probability of getting more methylated sites in certain genomic features ############


### Add 1 to all CpG loci to convert positions to 0 base to 1 base format ###
CpGsite_foot<-as.data.frame(y_foot) %>% 
  dplyr::select(Chr, Locus) %>% 
  mutate(Meth=0,
         Locus =Locus +1)

### Reformatting columns to include separate columns for chromosome and locus ###
DM_Trans_foot<- separate(DM_Trans_foot, sample, into = c("Chr", "Locus"), sep = "-")

### DM CpG subset associated with transplant site effect derived before, need to convert to 1 base format too ###
transDM_foot<-DM_Trans_foot%>% 
 mutate(Locus =as.numeric(Locus)+1)

#### Creating a dataframe with all the CpG sites, the ones that are differentially methylated with transplant site were "1" and nonmeth are "0" ###
df_trans <- CpGsite_foot %>%
  left_join(transDM_foot, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  dplyr::select(Chr,Locus, Meth)

  
### Adding genic features to the dataframe for each site ###

intron_subset_foot<-intron_subset_gene_id_df_foot %>% 
  mutate(feature = "intron") #25819 obs

exon_subset_foot<-exon_subset_gene_id_df_foot %>% 
  mutate(feature = "exon") #5114 obs

intergenic_subset_foot<-interg_subset_df_foot %>% 
  mutate(feature = "interg") #39339 obs

promoter_subset_foot<-promoter_subset_gene_id_df_foot %>% 
  mutate(feature = "promoter") #1729 obs

### Finding CpG subsets overlapped with exons ###  
exon<- CpGsite_foot %>%
  full_join(exon_subset_foot, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Finding CpG subsets overlapped with introns ###
intron <- CpGsite_foot %>%
  full_join(intron_subset_foot, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Finding CpG subsets overlapped with intergenic regions ###
intergenic <- CpGsite_foot %>%
  full_join(intergenic_subset_foot, by = c("Chr" = "seqnames", "Locus" = "start")) %>% 
    filter(!is.na(feature))

### Finding CpG subsets overlapped with promoter regions ###
promoter<-CpGsite_foot %>%
  full_join(promoter_subset_foot, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Combining all the CpG subsets ###
df_CpGsite_foot_all <- bind_rows(exon, intron, intergenic, promoter)

### Adding methylation information to the dataframe ###
df_combined_foot_trans <- df_trans%>% 
  full_join(df_CpGsite_foot_all , by =c("Chr","Locus")) %>% 
  select(-Meth.y) 

### Remove intergenic regions that overlap with promoter regions ###
  df_combined_foot_trans<-df_combined_foot_trans%>% 
    group_by(Chr, Locus) %>%
 filter(!(any(feature == "promoter") & feature == "interg")) %>%
  ungroup()

### Running GLM to test the effect of genomic feature on DM pattern (0 or 1) ###
glm_trans_foot <- glm(Meth.x ~ feature,
          data=df_combined_foot_trans, 
        family=binomial(link="logit"))
summary(glm_trans_foot)


### DM CpGs associated with origin site effect, separating columns into chromosome and locus ###
DM_Origin_foot<- separate(DM_Origin_foot, sample, into = c("Chr", "Locus"), sep = "-")

### Change loci to 0 base to 1 base ###
originDM_foot<-DM_Origin_foot%>% 
 mutate(Locus =as.numeric(Locus)+1)

### Adding methylation information for the CpGs, methylated as 1 and unmethylated as 0 ###
df_origin <- CpGsite_foot %>%
  left_join(originDM_foot, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  select(Chr,Locus, Meth)

df_combined_foot_origin <- df_origin %>% 
  full_join(df_CpGsite_foot_all , by =c("Chr","Locus")) %>% 
  select(-Meth.y) 

### Remove intergenic regions that overlap with promoter regions ###
df_combined_foot_origin <-df_combined_foot_origin%>% 
    group_by(Chr, Locus) %>%
 filter(!(any(feature == "promoter") & feature == "interg")) %>%
  ungroup()

### Running GLM to test the effect of feature on DM pattern (0 or 1) ###
glm_origin_foot <- glm(Meth.x ~ feature,
          data=df_combined_foot_origin, 
        family=binomial(link="logit"))
summary(glm_origin_foot)



############ GO Enrichment analysis ############

### Load in GO terms, downloaded from NCBI ###
go_terms<-
  read.delim(
    "/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/gene ontology/GCF_021869535.1_xbMytCali1.0.p_gene_ontology.gaf",
    header = FALSE,
    sep = "\t",
    # Specify tab delimiter
    comment.char = "!",
    # GAF files have comment lines starting with "!"
    stringsAsFactors = FALSE
  )  %>% # Avoid converting strings to factors) 
  dplyr::select(-V11, -V16, -V17)

### Changing column names ###
  colnames(go_terms) <- c("DB", "DB_Object_ID", "gene_id", "Qualifier", 
                         "GO_ID", "GO_Term", "Aspect", "DB_Object_Name", 
                         "Type", "DB_Object_Type", "Category","Taxon", 
                         "V13","V14") 

### Changing DM CpG to 1 base format ###
transDM_foot<-DM_Trans_foot %>% 
  mutate(Locus=as.numeric(Locus)+1)

originDM_foot<- DM_Origin_foot %>% 
  mutate(Locus=as.numeric(Locus)+1)

### Creating a dataframe with genomic feature & DM methylation classification, keep only entres with genes ###
df_trans_all_fix <- df_CpGsite_foot_all %>%
  left_join(transDM_foot, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  select(Chr,Locus, Meth,gene_id,feature) %>%  # Removing intergenic features that do not correspond with a gene
  na.omit()

### Extracting all DM CpGs ###
  DM_direc <- df_trans_all_fix %>%
  filter(Meth == 1) %>%
  na.omit()

### Extracting DM CpGs associated with both transplant and origin site effect ###
transplant_feature<-df_combined_foot_trans %>% filter(Meth.x ==1)%>%
  na.omit()

origin_feature<-df_combined_foot_origin %>% filter(Meth.x ==1)%>%
  na.omit()


DM_direc_origin<-left_join(DM_direc,
                           origin_feature,
                           by=c("gene_id","Locus","Chr"))

DM_direc_transplant<-left_join(DM_direc,
                            transplant_feature,
                            by=c("gene_id","Locus","Chr"))


### Matching genes with go terms ###
go_foot_trans_dm <- df_trans_all_fix %>% 
  left_join(go_terms, by = c("gene_id")) %>%  
  na.omit() 

### Count unique CpG sites per gene ###
cpg_per_gene<-go_foot_trans_dm  %>% 
  group_by(gene_id) %>% 
 summarize(unique_locus_count = n_distinct(paste(Locus, Chr))) %>% 
  na.omit()

### Extract unique gene IDs from GO terms ###
gene_id<-go_foot_trans_dm %>% select(gene_id) %>% 
  distinct()


### Filter out genes that have CpG coverage less than 3, save as cpg_sites_bias ###
cpg_sites_bias<- gene_id %>%
  left_join(cpg_per_gene, by = c("gene_id")) %>% 
  filter(unique_locus_count>=3) 

### GO Matching category file, which contains gene IDs and their corresponding GO IDs ###
category_file<-go_foot_trans_dm %>% 
  select(gene_id, GO_ID) %>% 
  distinct()%>%
  filter(gene_id %in% cpg_sites_bias$gene_id)%>% 
   distinct(gene_id, GO_ID)


############ For all GO TERMS extract the Parental level ############

### Loading libraries ###
library(GO.db)
library(AnnotationDbi) 

### Extract all the unique GO terms from category list ###
go_terms_cat <- unique(category_file$GO_ID)

### Extract ontology for each GO terms ###
go_info <- AnnotationDbi::select(GO.db, keys = go_terms_cat, columns = c("ONTOLOGY", "TERM"))

### Creating separate dataframe to store BP, MF, CC GO terms ###
bp<- go_info %>% 
  filter(ONTOLOGY == "BP") 
mf<-go_info %>% 
  filter(ONTOLOGY == "MF")
cc<-go_info %>% 
  filter(ONTOLOGY == "CC")

### BIOLOGICAL PROCESS ###
# Initialize a list to store parent terms
parent_terms <- list()

# Loop through each GO term and get the parent
for (go in bp$GOID) {
  # Attempt to retrieve the parent term (one level up)
  parent <- as.character(GO.db::GOBPPARENTS[[go]])
  
  # Debugging output
  if (is.null(parent) || length(parent) == 0) {
    cat("No parent term found for:", go, "\n")
    parent_terms[[go]] <- go  # Keep the original GO term if no parent is found
  } else {
    parent_terms[[go]] <- parent
  }
}

BP_parent_terms_df <- data.frame(
  GO_Term = rep(bp$GOID, sapply(parent_terms, length)),
  Parent_Term = unlist(parent_terms, use.names = FALSE),
  stringsAsFactors = FALSE) %>% 
  mutate(Ontology ="BP")


### CELLULAR COMPONENT ###
# Initialize a list to store parent terms
parent_terms <- list()

# Loop through each GO term and get the parent
for (go in cc$GOID) {
  # Attempt to retrieve the parent term
  parent <- as.character(GO.db::GOCCPARENTS[[go]])
  
  # Debugging output
  if (is.null(parent) || length(parent) == 0) {
    cat("No parent term found for:", go, "\n")
    parent_terms[[go]] <- go  # Keep the original GO term if no parent is found
  } else {
    parent_terms[[go]] <- parent
  }
}

CC_parent_terms_df <- data.frame(
  GO_Term = rep(cc$GOID, sapply(parent_terms, length)),
  Parent_Term = unlist(parent_terms, use.names = FALSE),
  stringsAsFactors = FALSE)%>% 
  mutate(Ontology ="CC")


  
### MOLECULAR FUNCTION ###
# Initialize a list to store parent terms
parent_terms <- list()

# Loop through each GO term and get the parent
for (go in mf$GOID) {
  # Attempt to retrieve the parent term
  parent <- as.character(GO.db::GOMFPARENTS[[go]])
  
  # Debugging output
  if (is.null(parent) || length(parent) == 0) {
    cat("No parent term found for:", go, "\n")
    parent_terms[[go]] <- go  # Keep the original GO term if no parent is found
  } else {
    parent_terms[[go]] <- parent
  }
}

MF_parent_terms_df <- data.frame(
  GO_Term = rep(mf$GOID, sapply(parent_terms, length)),
  Parent_Term = unlist(parent_terms, use.names = FALSE),
  stringsAsFactors = FALSE) %>% mutate(Ontology ="MF")


### Combine the BP, MF, CC data frames into one ###
all_parent_terms_df <- bind_rows(MF_parent_terms_df, CC_parent_terms_df, BP_parent_terms_df)

#### Rejoin with category file to get all the gene ID that matches those parental terms ###

category_file_parent<-category_file %>% 
  left_join(all_parent_terms_df, by=c("GO_ID" ="GO_Term"))

############################################################

### Only keep GO terms with at least 3 genes ###
gene_go_parent<-category_file_parent %>%
  group_by(Parent_Term) %>%
summarise(gene_count = n_distinct(gene_id)) %>% 
  filter(gene_count >= 3)

### Removing old GO terms, only including the parental GO terms ###
category_list_prep_parent <- gene_go_parent %>%
  left_join(category_file_parent, by = "Parent_Term") %>% 
  dplyr::select(-GO_ID)%>% 
   distinct(Parent_Term, gene_id, .keep_all = TRUE)

### Filtering out CpsG, include only those included in cpg_sites_bias (See above, only include genes with with at least 3 CpG coverage) ###
category_list <- cpg_sites_bias %>% 
  left_join(category_list_prep_parent, by = "gene_id") %>%  # Keep only genes in cpg_sites_bias
  na.omit() %>% 
  group_by(gene_id) %>%
  summarise(GO_ID = list(unique(Parent_Term)), .groups = 'drop') %>%  # Drop grouping afterwards
  deframe()

### The input bias correction file, which contains the gene id and the number of CpGs for that gene ###
cpg_sites_bias_parent<- category_list_prep_parent  %>% 
  left_join(cpg_sites_bias, by ="gene_id") %>% 
   dplyr::select(gene_id, unique_locus_count) %>% 
    distinct(gene_id, unique_locus_count) 

### Reformat bias correction file to vector format ###
cpg_vector <- setNames(cpg_sites_bias_parent$unique_locus_count,
                          cpg_sites_bias_parent$gene_id)

### Create dataframe with gene id, CpG count for each locus, methylation information, and GO annotations ###
go_foot_trans_dm_parent <- cpg_sites_bias_parent %>% 
  left_join(go_foot_trans_dm, by ="gene_id")

  
### Separating unmethylated and methylated CpGs ###
DMG<-go_foot_trans_dm_parent%>% 
  filter(Meth ==1) %>% 
   dplyr::select(gene_id, GO_ID) %>% 
  distinct()

UD<-go_foot_trans_dm_parent %>% 
  filter(Meth ==0) %>% 
   dplyr::select(gene_id, GO_ID) %>% 
  distinct()

### extract DM and Unmeth gene ID ###
dm_genes <- DMG$gene_id
ud_genes <- UD$gene_id

### Combine all genes into a single vector ###
all_genes <- union(dm_genes, ud_genes)

### Create a named vector with DM genes as 1 and non-DM genes as 0 ###
gene_vector <- setNames(
  as.integer(all_genes %in% de_genes),
  all_genes
)
 
### Running GO enrichment analysis with bias data ###
pwf <- nullp(gene_vector, bias.data=cpg_vector, plot.fit =FALSE)
head(pwf)

plotPWF(pwf = pwf, binsize = 200) 

GO.wall <- goseq(pwf, gene2cat = category_list, method = "Wallenius", use_genes_without_cat = FALSE)

overrep_go<-GO.wall %>% filter(over_represented_pvalue<0.05)
enriched_go_ids <- overrep_go$category

 ### FDR corrected ###
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH" ) < 0.05]
 head(enriched.GO)
 #character(0)

