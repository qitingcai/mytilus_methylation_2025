### gill specific ###
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

### Set working directory, include all (top5 for each treatments) gill coverage files ###
setwd("/Users/qcai/Documents/UCSC/Kelley_Lab/mytilus/cg_coverage_files/tissue_specific/gill")

### Load metadata file that include sample names, load all coverage files into R ###
meta_data_gill<-read.delim("gill_metadata.txt", row.names = "sample", stringsAsFactors = FALSE) 
Sample_gill <- row.names(meta_data_gill)
files_gill <- paste0(Sample_gill,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

### EdgeR function readBismark2DGE reads all the files and collates the counts for all the sample into one data object ###
yall <- readBismark2DGE(files_gill, sample.names=Sample_gill)

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
keep_gill <- rowSums(Coverage >= n) >= 14
table(keep_foot)

### DGEList object is subsetted to retain only the filtered loci ###
y_gill <- yall[keep_gill,, keep.lib.sizes=FALSE]

### Normalization - set the library sizes for each sample to be the average of the total read counts for the methylated and unmethylated libraries ###
TotalLibSize <- y_gill$samples$lib.size[Methylation=="Me"] +
  +                 y_gill$samples$lib.size[Methylation=="Un"]
y_gill$samples$lib.size <- rep(TotalLibSize, each=2)
 y_gill$samples

### Compute the corresponding methylation summary from the methylated and unmethylated counts ###
Me_gill <- y_gill$counts[, Methylation=="Me"]
Un_gill <- y_gill$counts[, Methylation=="Un"]

### Calculating a methylation proportion matrix ###
prop_meth_matrix_gill <- Me_gill/(Me_gill+Un_gill)

### We want to use the exposed site as reference level for both transplant and origin site effect analyses, need to relevel for transplant site effect ###
meta_data_gill$transplant_site<- relevel(factor(meta_data_gill$transplant_site), ref = "exposed")

### Create a design matrix, using origin site, transplant site, and final shell length as fixed effects ###
designSL_gill <- model.matrix(~0+ origin_site + transplant_site +
                           LENGTH_FINAL..mm., 
                         data=meta_data_gill)

### Expand to the full design matrix modeling the sample and methylation effects ###
design_gill <- modelMatrixMeth(designSL_gill)

### Dispersion estimation ###
y_gill <- estimateDisp(y_gill, design = design_gill, robust = TRUE)

### Create the BCV plot ### 
plotBCV(y_gill)


### Testing for differentially methylated CpG loci ###

### fit NB GLMs for all the CpG loci using the glmFit function in edgeR ###
fit_gill <- glmFit(y_gill, design_gill)

### Testing for differentially methylated CpG sites between different treatment groups using likelihood ratio tests ###
### Origin site effects ###
 contr_origin_gill <- makeContrasts(Origin = origin_siteprotected-origin_siteexposed, levels = design_gill)
lrt_origin_gill <- glmLRT(fit_gill, contrast=contr_origin_gill)
summary( decideTests(lrt_origin_gill) )

### Transplant site effects (exposed is the reference level) ###
contr_trans_gill <- makeContrasts(Transplant = transplant_siteprotected, levels = design_gill)
lrt_trans_gill <- glmLRT(fit_gill, contrast=contr_trans_gill)
summary( decideTests(lrt_trans_gill) )


### Wrangle data for generating volano plot of CpG DM ###
### Correct p-values using BH method ###
lrt_origin_gill$table$FDR <- p.adjust( lrt_origin_gill$table$PValue, method = "BH" )
lrt_trans_gill$table$FDR <- p.adjust( lrt_trans_gill$table$PValue, method = "BH" )

# Apply logical variable for significant DM ###
lrt_origin_gill$table$Sig <- ifelse(lrt_origin_gill$table$FDR< 0.05, TRUE, FALSE )
lrt_trans_gill$table$Sig <- ifelse(lrt_trans_gill$table$FDR< 0.05, TRUE, FALSE )

# Apply binary variable for hyper/hypomethylation or "Up" vs "Down" ###
lrt_origin_gill$table$Dir <- ifelse(lrt_origin_gill$table$logFC > 0, "Up", "Down" )
lrt_trans_gill$table$Dir <- ifelse( lrt_trans_gill$table$logFC > 0, "Up", "Down" )

# Create combined term for significance and fold-change direction of diff meth ###
lrt_origin_gill$table$Sig_Dir <- paste( lrt_origin_gill$table$Sig,
                             lrt_origin_gill$table$Dir,
                             sep = "_" )
lrt_trans_gill$table$Sig_Dir <- paste( lrt_trans_gill$table$Sig,
                             lrt_trans_gill$table$Dir,
                             sep = "_" )

### Create CpG ID and treatment variables ###
lrt_origin_gill$table$Treat <- "Origin"
lrt_trans_gill$table$Treat <- "Transplant"

lrt_origin_gill$table$sample <- rownames(lrt_origin_gill$table)
lrt_trans_gill$table$sample <- rownames(lrt_trans_gill$table)

### Merge origin site and transplant site coefficients ###
all_CpG_dm <- rbind( lrt_origin_gill$table,
                     lrt_trans_gill$table )

### Saving DM sites for transplant site effect and origin site effect as new objects ###
DM_Trans_gill<-all_CpG_dm %>% filter(Sig =="TRUE",Treat =="Transplant")
DM_Origin_gill<-all_CpG_dm %>% filter(Sig =="TRUE", Treat =="Origin")

### Origin site and transplant site associated DM CpG volcano plots ###
ggplot( data = all_CpG_dm, aes(y = -log( as.numeric( FDR ) ), x = as.numeric( logFC ),
                               color = Sig_Dir ) ) +
    geom_point() +
    theme_classic( base_size = 40 ) +
    theme( legend.position = "none",
           strip.background = element_blank() ) +
    scale_color_manual( values = c( "black", "black", "blue", "red" ) ) +
    facet_grid( Treat ~ . , scale = "free") +
    labs( x = "Gill Diff meth", y = "-log FDR" )


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
  slice_min(order_by = start, n = 1, with_ties =FALSE) #Only return one exon when there are overlapping first exons

### For negative strands, extract the first exons using the largest start position with slice_max, (going from right to left) ###
first_exon_neg <- as.data.frame(exon_neg_strand) %>% 
  group_by(gene) %>% 
  slice_max(order_by = start, n = 1, with_ties =FALSE) #

### Combine both positive and negative strand first exons ###
first_exon_GRange <- bind_rows(first_exon_pos, first_exon_neg)

### Remove potential duplicated lines ###
first_exon_GRange <- first_exon_GRange %>% distinct()

### Create new GRange object with only first exons ###
first_exon_GRange <- GRanges(
  seqnames = Rle(first_exon_GRange$seqnames),
  ranges = IRanges(start = first_exon_GRange$start, end = first_exon_GRange$end),
  strand = first_exon_GRange$strand,
  gene_id =first_exon_GRange$gene
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
  1000,  # If positive strand and start is 1, max end at 1000, only for cases where first exon starts at position <1000
  ifelse(
    strand(first_exon_GRange) == "+",
    start(first_exon_GRange) - 1,  # For positive strand, end before the first position of first exon
    end(first_exon_GRange)+1001  # For negative strand, end 1000bp downstream of the end of first exon, going from left to right
  )
)

### Create the GRanges object for promoter regions ###
promoters_Grange <- GRanges(
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
prop_meth_matrix_gill <- as.data.frame(prop_meth_matrix_gill)

### Reformating ###
prop_meth_matrix_gill <- prop_meth_matrix_gill %>%
 rownames_to_column(var = "name") #make the rows become first column

CpGs_gill <- prop_meth_matrix_gill  %>%
  separate(col = 1, into = c("seqid", "start"), sep = "-") #reformatting the CpG location identifiers

CpGs_gill$start<-as.numeric(CpGs_gill$start)

CpGs_gill<-CpGs_gill %>% 
  mutate(start = start + 1,
         end = start +1) #changing from 0 to 1 base for cpg positions, because the annotation files are 1 base, but the output from bismark is 0 base

### Reorder columns to place 'end' in the third column for easier viewing ###
CpGs_gill <- CpGs_gill[, c(1, 2, ncol(CpGs_gill), 3:(ncol(CpGs_gill)-1))]

### Make cpg a GRange object for all CpG sites ###
CpGs_gill<-makeGRangesFromDataFrame(CpGs_gill)

### Subsetting CpG GRange objects according to genic features ###
### Exon subset ###
exon_subset_gill<-subsetByOverlaps(CpGs_gill, exon_GRange)

### Adding gene id to each entry ###
exon_gff <- subset(all_GRange, type == "exon") #need to be GRange object to subset
mcols(exon_gff)
exon_id_overlaps_gill <- findOverlaps(exon_subset_gill, exon_gff)

### Extract gene_id for matching exons from the GFF file ###
gene_ids_gill <- mcols(exon_gff)$gene[subjectHits(exon_id_overlaps_gill)]
query_hits_gill <- queryHits(exon_id_overlaps_gill)
# Add gene_id only to those rows in exon_subset that have corresponding overlaps

### Initialize with NA for those rows that have no gene overlap ###
exon_subset_gene_id_gill <- rep(NA, length(exon_subset_gill))
# Assign the gene_ids where there is a match
exon_subset_gene_id_gill[query_hits_gill] <- gene_ids_gill

### Add gene_id as a metadata column to exon_subset df ###
mcols(exon_subset_gill)$gene_id <- exon_subset_gene_id_gill
exon_subset_gene_id_df_gill<-as.data.frame(exon_subset_gill)
exon_subset_gene_id_df_gill$combined <- paste(exon_subset_gene_id_df_gill$seqnames, exon_subset_gene_id_df_gill$start, sep = "-") #reformat


### Intron subset ###
intron_subset_gill<-subsetByOverlaps(CpGs_gill, intron_GRange)
#adding gene id
intron_gff <- subset(all_GRange, type == "intron") #need to be GRange object to subset
mcols(intron_gff)
intron_id_overlaps_gill <- findOverlaps(intron_subset_gill, intron_gff)

### Extract gene_id for matching exons from the GFF file ###
gene_ids_gill<- mcols(intron_gff)$gene[subjectHits(intron_id_overlaps_gill)]
query_hits <- queryHits(intron_id_overlaps_gill)


### Add gene_id only to those rows in exon_subset that have corresponding overlaps ###
### Initialize with NA for those rows that have no overlap ###
intron_subset_gene_id_gill <- rep(NA, length(intron_subset_gill))

### Assign the gene_ids where there is a match ###
intron_subset_gene_id_gill[query_hits] <- gene_ids_gill

### Add gene_id as a metadata column to intron_subset df ###
mcols(intron_subset_gill)$gene_id <- intron_subset_gene_id_gill
intron_subset_gene_id_df_gill<-as.data.frame(intron_subset_gill)
intron_subset_gene_id_df_gill$combined <- paste(intron_subset_gene_id_df_gill$seqnames, intron_subset_gene_id_df_gill$start, sep = "-") #reformat

### Find promoter overlaps with CpG using the adjusted ranges ###
promoter_GRange_gill <- findOverlaps(CpGs_gill, promoters_Grange)
promoter_subset_gill<-subsetByOverlaps(CpGs_gill, promoters_Grange)


### Adding gene id ###
promoter_gff <- promoters_Grange #need to be GRange object to subset
promoter_id_overlaps_gill <- findOverlaps(promoter_subset_gill, promoter_gff)

### Extract gene_id for matching promoter region (genes associated with first exons, see above, from the GFF file ###
gene_ids_gill<- mcols(promoter_gff)$gene[subjectHits(promoter_id_overlaps_gill)]
query_hits <- queryHits(promoter_id_overlaps_gill)

### Add gene_id only to those rows in promoter_subset_foot that have corresponding overlaps ###
# Initialize with NA for those rows that have no overlap
promoter_subset_gene_id_gill <- rep(NA, length(promoter_subset_gill))

### Assign the gene_ids where there is a match ###
promoter_subset_gene_id_gill[query_hits] <- gene_ids_gill

### Add gene_id as a metadata column to promoter subset df ###
mcols(promoter_subset_gill)$gene_id <- promoter_subset_gene_id_gill
promoter_subset_gene_id_df_gill<-as.data.frame(promoter_subset_gill)
promoter_subset_gene_id_df_gill$combined <- paste(promoter_subset_gene_id_df_gill$seqnames, promoter_subset_gene_id_df_gill$start, sep = "-") #reformat


### FINDING INTERGENIC REGIONS ###
### Filter for all genic regions (exon, intron, UTRs) ###
non_interg_GRange <- subset(all_GRange, type %in% c("exon", "intron", "five_prime_UTR", "three_prime_UTR"))

### Find overlaps between CpG and non-intergenic regions (genic regions) ###
noninterg_overlap <- findOverlaps(CpGs_gill, non_interg_GRange)

### Get the indices of the CpG regions that overlap with genic regions ###
overlapping_indices <- queryHits(noninterg_overlap)

### Remove these overlapping CpG regions to get intergenic regions ###
interg_subset <- CpGs_gill[-overlapping_indices]

### Convert to data frame and make a Grange Object, genes are NAs ###
interg_subset_df_gill <- as.data.frame(interg_subset) %>% mutate(gene_id =NA)
interg_subset_df_gill$combined <- paste(interg_subset_df_gill$seqnames, interg_subset_df_gill$start, sep = "-")


############ Running GLM to assess the probability of getting more methylated sites in certain genomic features ############


### Add 1 to all CpG loci to convert positions to 0 base to 1 base format ###

CpGsite_gill<-as.data.frame(y_gill) %>% 
  dplyr::select(Chr, Locus) %>% 
  mutate(Meth=0,
         Locus =Locus +1)

### Reformatting columns to include separate columns for chromosome and locus ###
DM_Trans_gill<- separate(DM_Trans_gill, sample, into = c("Chr", "Locus"), sep = "-")

### DM CpG subset associated with transplant site effect derived before, need to convert to 1 base format too ###
transDM_gill<-DM_Trans_gill%>% 
 mutate(Locus =as.numeric(Locus)+1)

#### Creating a dataframe with all the CpG sites, the ones that are differentially methylated with transplant site were "1" and nonmeth are "0" ###
df_trans <- CpGsite_gill %>%
  left_join(transDM_gill, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  dplyr::select(Chr,Locus, Meth)

### Adding genic features to the dataframe for each site ###

intron_subset_gill<-intron_subset_gene_id_df_gill %>% 
  mutate(feature = "intron") #28378 obs

exon_subset_gill<-exon_subset_gene_id_df_gill %>% 
  mutate(feature = "exon") #5645 obs

intergenic_subset_gill<-interg_subset_df_gill %>% 
  mutate(feature = "interg") #43382 obs

promoter_subset_gill<-promoter_subset_gene_id_df_gill %>% 
  mutate(feature = "promoter") #1990 obs

### Finding CpG subsets overlapped with exons ###  
exon<- CpGsite_gill %>%
  full_join(exon_subset_gill, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Finding CpG subsets overlapped with introns ###
intron <- CpGsite_gill %>%
  full_join(intron_subset_gill, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Finding CpG subsets overlapped with intergenic regions ###
intergenic <- CpGsite_gill %>%
  full_join(intergenic_subset_gill, by = c("Chr" = "seqnames", "Locus" = "start")) %>% 
    filter(!is.na(feature))

### Finding CpG subsets overlapped with promoter regions ###
promoter<-CpGsite_gill %>%
  full_join(promoter_subset_gill, by = c("Chr" = "seqnames", "Locus" = "start")) %>%
na.omit()

### Combining all the CpG subsets ###
df_CpGsite_gill_all <- bind_rows(exon, intron, intergenic, promoter)

### Adding methylation information to the dataframe ###
df_combined_gill_transplant <- df_trans%>% 
  full_join(df_CpGsite_gill_all , by =c("Chr","Locus")) %>% 
  dplyr::select(-Meth.y) 


### Running GLM to test the effect of feature on DM pattern (0 or 1) ###
glm_transplant_gill <- glm(Meth.x ~ feature,
          df_combined_gill_transplant , 
          family=binomial(link="logit"))
summary(glm_transplant_gill)


### DM CpGs associated with origin site effect, separating columns into chromosome and locus ###
DM_Origin_gill<- separate(DM_Origin_gill, sample, into = c("Chr", "Locus"), sep = "-")

### Change loci to 0 base to 1 base ###
originDM_gill<-DM_Origin_gill%>% 
 mutate(Locus =as.numeric(Locus)+1)

### Adding methylation information for the CpGs, methylated as 1 and unmethylated as 0 ###
df_origin <- CpGsite_gill %>%
  left_join(originDM_gill, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  dplyr::select(Chr,Locus, Meth)

df_combined_gill_origin <- df_origin%>% 
  full_join(df_CpGsite_gill_all , by =c("Chr","Locus")) %>% 
  dplyr::select(-Meth.y) 

### Running GLM to test the effect of genomic feature on DM pattern (0 or 1) ###
glm_origin_gill <- glm(Meth.x ~ feature,
          df_combined_gill_origin, 
          family=binomial(link="logit"))
summary(glm_origin_gill)



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
transDM_gill<- DM_Trans_gill %>% 
  mutate(Locus = as.numeric(Locus)+1) 

originDM_gill<- DM_Origin_gill %>% 
  mutate(Locus = as.numeric(Locus)+1) 

### Creating a dataframe with genomic feature & DM methylation classification, keep only entres with genes ###
df_trans_all_fix <- df_CpGsite_gill_all %>%
  left_join(transDM_gill, by = c("Chr", "Locus")) %>%
  mutate(Meth = if_else(!is.na(logFC), 1, Meth)) %>%
  select(Chr,Locus, Meth,gene_id) %>%  # Removing intergenic features that do not correspond with a gene
  na.omit()

### Matching genes with go terms ###
go_gill_trans_dm <- df_trans_all_fix %>% 
  left_join(go_terms, by = c("gene_id")) %>%  
  na.omit() 

### Count unique CpG sites per gene ###
cpg_per_gene<-go_gill_trans_dm  %>% 
  group_by(gene_id) %>% 
 summarize(unique_locus_count = n_distinct(paste(Locus, Chr))) %>% 
  na.omit()

### Extract unique gene IDs from GO terms ###
gene_id<-go_gill_trans_dm %>% select(gene_id) %>% 
  distinct()

### Filter out genes that have CpG coverage less than 3, save as cpg_sites_bias ###
cpg_sites_bias<- gene_id %>%
  left_join(cpg_per_gene, by = c("gene_id")) %>% 
  filter(unique_locus_count>=3) 

### GO Matching category file, which contains gene IDs and their corresponding GO IDs ###
category_file<-go_gill_trans_dm %>% 
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
go_info <- AnnotationDbi::select(GO.db, keys = go_terms_cat, columns = c("ONTOLOGY", "TERM"))

### Extract ontology for each GO terms ###
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
  # Attempt to retrieve the parent term
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
  stringsAsFactors = FALSE) %>%  mutate(Ontology ="BP")


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
  stringsAsFactors = FALSE
) %>% 
  mutate(Ontology ="MF")


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
  left_join(category_file_parent, by = "Parent_Term") #go terms

### Filtering out CpsG, include only those included in cpg_sites_bias (See above, only include genes with with at least 3 CpG coverage) ###
category_list <- category_list_prep_parent %>%  
  group_by(gene_id) %>%
  summarise(GO_ID = list(unique(Parent_Term)), .groups = 'drop') %>%  # Drop grouping afterwards
  deframe()

### The input bias correction file, which contains the gene id and the number of CpGs for that gene ### 
cpg_sites_bias_parent<- category_list_prep_parent %>% 
  left_join(cpg_sites_bias, by ="gene_id") %>% 
  dplyr::select(gene_id, unique_locus_count) %>% 
    distinct(gene_id, unique_locus_count) 


### Reformat bias correction file to vector format ###

cpg_vector <- setNames(cpg_sites_bias_parent$unique_locus_count,
                          cpg_sites_bias_parent$gene_id)


### Create dataframe with gene id, CpG count for each locus, methylation information, and GO annotations ###
go_gill_trans_dm_parent <- cpg_sites_bias_parent %>% 
  left_join(go_gill_trans_dm, by ="gene_id")


### Separating unmethylated and methylated CpGs ###
DMG<-go_gill_trans_dm_parent %>% 
  filter(Meth ==1) %>% 
  dplyr::select(gene_id, GO_ID) %>% 
  distinct()

UD<-go_gill_trans_dm_parent %>% 
  filter(Meth ==0) %>% 
   dplyr::select(gene_id, GO_ID) %>% 
  distinct()

### extract DM and Unmeth gene ID ###
de_genes <- DMG$gene_id
ud_genes <- UD$gene_id


### Combine all genes into a single vector ###
all_genes <- union(de_genes, ud_genes)

### Create a named vector with DM genes as 1 and non-DM genes as 0 ###
gene_vector <- setNames(
  as.integer(all_genes %in% de_genes),
  all_genes
)

pwf <- nullp(gene_vector, bias.data=cpg_vector, plot.fit =FALSE)
head(pwf)

plotPWF(pwf = pwf, binsize = 200) 

GO.wall <- goseq(pwf, gene2cat = category_list, method = "Wallenius", use_genes_without_cat = FALSE)

overrep_go<-GO.wall %>% filter(over_represented_pvalue<0.05)
 enriched_go_ids <- overrep_go$category

 
 ### FDR corrected ###
 enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH" ) < .05]
 head(enriched.GO) #character(0)