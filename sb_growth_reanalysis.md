sb_growth_reanalysis
================
Sam Bogan
2026-02-26

This is a reanalysis and adaptation of Tina’s code by Sam in response to
the comment by reviewer 2 that we should evaluate initial shell size
and/or growth rate as a covariate in our DM tests.

The first thing we’ll want to do is see if the association between
initial shell size and growth rate is flat (linear growth), positive
(exponential growth) or non-linear/negative (logarithmic growth).

## Look at shell length and growth

``` r
# Read in shell metadata
shell_df <- read.csv("Data/full_sample_sbcorr.csv")

# Calc linear growth rate and make one row for each mussel individual
shell_df <- shell_df %>%
  dplyr::mutate(
    DATE_INIT = mdy(DATE_INIT),
    DATE_SAMPLED = mdy(DATE_SAMPLED)
    ) %>%
  mutate(
    days_elapsed = as.numeric(DATE_SAMPLED - DATE_INIT),
    mm_per_d = (LENGTH_FINAL..mm. - LENGTH_INIT..mm.) / days_elapsed,
    prop_mm_per_d = ((LENGTH_FINAL..mm. - LENGTH_INIT..mm.) / days_elapsed) / LENGTH_INIT..mm.
  ) %>%
  distinct(sample.ID, .keep_all = TRUE)

# Plot initial length vs linear growth rate
# growth is linearbut doesn't correlate well with initial shell size
ggplot(data = shell_df,
       aes(x = `LENGTH_INIT..mm.`, y = mm_per_d)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(x = "Initial shell length (mm)", y = "Linear growth rate (mm per day)")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Now proportional growth
ggplot(data = shell_df,
       aes(x = `LENGTH_INIT..mm.`, y = prop_mm_per_d)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(x = "Initial shell length (mm)", y = "Proportional growth per day")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
# What about final shell length and growth rate?
# Slightly better correlation, and expected positive relationship...
# Higher growth rate results in bigger final shell length
ggplot(data = shell_df,
       aes(y = `LENGTH_FINAL..mm.`, x = mm_per_d)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(y = "Final shell length (mm)", x = "Linear growth rate (mm per day)")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
# Now proportional growth
ggplot(data = shell_df,
       aes(y = `LENGTH_FINAL..mm.`, x = prop_mm_per_d)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(y = "Final shell length (mm)", x = "Proportional growth per day")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

``` r
# What about correlation between initial and final size?
ggplot(data = shell_df,
       aes(x = `LENGTH_INIT..mm.`, y = LENGTH_FINAL..mm.)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(x = "Initial shell length (mm)", y = "Final shell length (mm per day)")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
# No correlation is tight enough to say that final shell length is a proxy for initial

# Proportional growth rate is non-independent of initial length, so I'm going to avoid using that

# Three models make sense and do not run into non-independence issues
# 1. initial length (basically the same as final growth according to correlation)
# 2. growth rate
# 3. initial length + mm growth rate

# There is no need to fit a model with proportional growth bc 
# the slope of mm growth across initial size = 0
```

Runs of Tina’s DM code with the three models above (one old; two new)

``` r
# Load DM packages
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 4.2.3

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-6.1

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(SYNCSA)
library(mice)
```

    ## 
    ## Attaching package: 'mice'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

``` r
library(ape)
```

    ## Warning: package 'ape' was built under R version 4.2.3

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

``` r
library(rtracklayer)
```

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:mice':
    ## 
    ##     cbind, rbind

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: GenomeInfoDb

``` r
library(genomation)
```

    ## Loading required package: grid

    ## Warning: replacing previous import 'Biostrings::pattern' by 'grid::pattern'
    ## when loading 'genomation'

``` r
library(plyranges)
```

    ## 
    ## Attaching package: 'plyranges'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, n, n_distinct

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(GenomicRanges)
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 4.2.3

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.2.3

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(ggthemes)
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     desc

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(coin)
```

    ## Loading required package: survival

    ## 
    ## Attaching package: 'coin'

    ## The following objects are masked from 'package:rstatix':
    ## 
    ##     chisq_test, friedman_test, kruskal_test, sign_test, wilcox_test

``` r
library(tidyverse)
```

## Start with gill

``` r
### Load metadata file that include sample names, load all coverage files into R ###
setwd("Data/gill_coverage_files/")
meta_data_gill<-read.delim("gill_metadata.txt", row.names = "sample", stringsAsFactors = FALSE)
Sample_gill <- row.names(meta_data_gill)
files_gill <- paste0("gill_samples/", Sample_gill,".CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov")

### EdgeR function readBismark2DGE reads all the files and collates the counts for all the sample into one data object ###
yall <- readBismark2DGE(files_gill, sample.names=Sample_gill)
```

    ## Reading gill_samples/14B-G_S25.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/15G-G_S24.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/18W-G_S41.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/21B-G_S7.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/29W-G_S76.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/29Y-G_S68.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/2B-G_S17.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/40Y-G_S63.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/45G-G_S69.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/47W-G_S80.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/50B-G_S83.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/51B-G_S18.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/55B-G_S40.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/55Y-G_S10.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/56Y-G_S53.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/67B-G_S8.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/67W-G_S86.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/75Y-G_S9.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/76Y-G_S85.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Reading gill_samples/90W-G_S67.CpG_report.merged_CpG_evidence.cov.CpG_report.merged_CpG_evidence.cov 
    ## Hashing ...
    ## Collating counts ...

``` r
### Check dimension and the count matrix ###
dim(yall)
```

    ## [1] 2485629      40

``` r
head(yall$counts)
```

    ##                      14B-G_S25-Me 14B-G_S25-Un 15G-G_S24-Me 15G-G_S24-Un
    ## NW_026262581.1-97825            2            5            0            4
    ## NW_026262581.1-97832            0            7            0            4
    ## NW_026262581.1-97836            0            7            0            4
    ## NW_026262581.1-97843            0            7            0            4
    ## NW_026262581.1-97856            0            7            0            4
    ## NW_026262581.1-97859            0            7            0            4
    ##                      18W-G_S41-Me 18W-G_S41-Un 21B-G_S7-Me 21B-G_S7-Un
    ## NW_026262581.1-97825            0            0           0          13
    ## NW_026262581.1-97832            0            0           0          13
    ## NW_026262581.1-97836            0            0           0          13
    ## NW_026262581.1-97843            0            0           0          13
    ## NW_026262581.1-97856            0            0           0          13
    ## NW_026262581.1-97859            0            0           0          13
    ##                      29W-G_S76-Me 29W-G_S76-Un 29Y-G_S68-Me 29Y-G_S68-Un
    ## NW_026262581.1-97825            1            3            0            2
    ## NW_026262581.1-97832            0            4            0            2
    ## NW_026262581.1-97836            0            4            0            2
    ## NW_026262581.1-97843            0            5            0            2
    ## NW_026262581.1-97856            0            5            0            2
    ## NW_026262581.1-97859            0            5            0            2
    ##                      2B-G_S17-Me 2B-G_S17-Un 40Y-G_S63-Me 40Y-G_S63-Un
    ## NW_026262581.1-97825           0           4            0            9
    ## NW_026262581.1-97832           0           4            0            9
    ## NW_026262581.1-97836           0           4            0            9
    ## NW_026262581.1-97843           0           4            0            9
    ## NW_026262581.1-97856           0           4            0            9
    ## NW_026262581.1-97859           0           4            0            9
    ##                      45G-G_S69-Me 45G-G_S69-Un 47W-G_S80-Me 47W-G_S80-Un
    ## NW_026262581.1-97825            0            4            0            4
    ## NW_026262581.1-97832            0            4            0            4
    ## NW_026262581.1-97836            0            4            0            4
    ## NW_026262581.1-97843            0            4            0            4
    ## NW_026262581.1-97856            0            4            0            4
    ## NW_026262581.1-97859            0            4            0            4
    ##                      50B-G_S83-Me 50B-G_S83-Un 51B-G_S18-Me 51B-G_S18-Un
    ## NW_026262581.1-97825            0           10            0            4
    ## NW_026262581.1-97832            0           10            0            4
    ## NW_026262581.1-97836            0           10            0            4
    ## NW_026262581.1-97843            0           10            0            4
    ## NW_026262581.1-97856            0           10            0            4
    ## NW_026262581.1-97859            0           10            0            4
    ##                      55B-G_S40-Me 55B-G_S40-Un 55Y-G_S10-Me 55Y-G_S10-Un
    ## NW_026262581.1-97825            1            7            0            2
    ## NW_026262581.1-97832            0            8            0            2
    ## NW_026262581.1-97836            0            8            0            2
    ## NW_026262581.1-97843            0            8            0            2
    ## NW_026262581.1-97856            0            8            0            2
    ## NW_026262581.1-97859            0            8            0            2
    ##                      56Y-G_S53-Me 56Y-G_S53-Un 67B-G_S8-Me 67B-G_S8-Un
    ## NW_026262581.1-97825            0            2           1          11
    ## NW_026262581.1-97832            0            2           0          12
    ## NW_026262581.1-97836            0            2           0          12
    ## NW_026262581.1-97843            0            2           0          12
    ## NW_026262581.1-97856            0            2           0          12
    ## NW_026262581.1-97859            0            2           0          12
    ##                      67W-G_S86-Me 67W-G_S86-Un 75Y-G_S9-Me 75Y-G_S9-Un
    ## NW_026262581.1-97825            0            3           0           2
    ## NW_026262581.1-97832            0            3           0           2
    ## NW_026262581.1-97836            0            3           0           2
    ## NW_026262581.1-97843            0            3           0           2
    ## NW_026262581.1-97856            0            3           0           2
    ## NW_026262581.1-97859            0            3           0           2
    ##                      76Y-G_S85-Me 76Y-G_S85-Un 90W-G_S67-Me 90W-G_S67-Un
    ## NW_026262581.1-97825            0           14            0            4
    ## NW_026262581.1-97832            0           14            0            4
    ## NW_026262581.1-97836            0           14            0            4
    ## NW_026262581.1-97843            0           14            0            4
    ## NW_026262581.1-97856            0           14            0            4
    ## NW_026262581.1-97859            0           14            0            4

``` r
### Save a data frame for downstream analyses ###
yall_df<-as.data.frame(yall)

### Sum up the counts of methylated and unmethylated reads to get the total read coverage at each CpG site for each sample ###
Methylation <- gl(2, 1, ncol(yall), labels = c("Me", "Un"))
Me <- yall$counts[ , Methylation == "Me" ]
Un <- yall$counts[ , Methylation == "Un" ]
Coverage <- Me + Un
# Prefiltered total # of CpGs sequenced: colSums(Coverage > 0, na.rm = TRUE) 
head(Coverage)
```

    ##                      14B-G_S25-Me 15G-G_S24-Me 18W-G_S41-Me 21B-G_S7-Me
    ## NW_026262581.1-97825            7            4            0          13
    ## NW_026262581.1-97832            7            4            0          13
    ## NW_026262581.1-97836            7            4            0          13
    ## NW_026262581.1-97843            7            4            0          13
    ## NW_026262581.1-97856            7            4            0          13
    ## NW_026262581.1-97859            7            4            0          13
    ##                      29W-G_S76-Me 29Y-G_S68-Me 2B-G_S17-Me 40Y-G_S63-Me
    ## NW_026262581.1-97825            4            2           4            9
    ## NW_026262581.1-97832            4            2           4            9
    ## NW_026262581.1-97836            4            2           4            9
    ## NW_026262581.1-97843            5            2           4            9
    ## NW_026262581.1-97856            5            2           4            9
    ## NW_026262581.1-97859            5            2           4            9
    ##                      45G-G_S69-Me 47W-G_S80-Me 50B-G_S83-Me 51B-G_S18-Me
    ## NW_026262581.1-97825            4            4           10            4
    ## NW_026262581.1-97832            4            4           10            4
    ## NW_026262581.1-97836            4            4           10            4
    ## NW_026262581.1-97843            4            4           10            4
    ## NW_026262581.1-97856            4            4           10            4
    ## NW_026262581.1-97859            4            4           10            4
    ##                      55B-G_S40-Me 55Y-G_S10-Me 56Y-G_S53-Me 67B-G_S8-Me
    ## NW_026262581.1-97825            8            2            2          12
    ## NW_026262581.1-97832            8            2            2          12
    ## NW_026262581.1-97836            8            2            2          12
    ## NW_026262581.1-97843            8            2            2          12
    ## NW_026262581.1-97856            8            2            2          12
    ## NW_026262581.1-97859            8            2            2          12
    ##                      67W-G_S86-Me 75Y-G_S9-Me 76Y-G_S85-Me 90W-G_S67-Me
    ## NW_026262581.1-97825            3           2           14            4
    ## NW_026262581.1-97832            3           2           14            4
    ## NW_026262581.1-97836            3           2           14            4
    ## NW_026262581.1-97843            3           2           14            4
    ## NW_026262581.1-97856            3           2           14            4
    ## NW_026262581.1-97859            3           2           14            4

``` r
## Calculating # of unique CpGs 
# covered_per_cpg <- rowSums(Coverage > 0, na.rm = TRUE)
# sum(covered_per_cpg > 0) # 1946085

### Filtering to only include samples with a coverage of at least 3, for 66% of samples, total of 14 samples ###
n=3
keep_gill <- rowSums(Coverage >= n) >= 14
table(keep_gill)
```

    ## keep_gill
    ##   FALSE    TRUE 
    ## 2408993   76636

``` r
### DGEList object is subsetted to retain only the filtered loci ###
y_gill <- yall[keep_gill,, keep.lib.sizes=FALSE]

### Normalization - set the library sizes for each sample to be the average of the total read counts for the methylated and unmethylated libraries ###
TotalLibSize <- y_gill$samples$lib.size[Methylation=="Me"] +
  +                 y_gill$samples$lib.size[Methylation=="Un"]
y_gill$samples$lib.size <- rep(TotalLibSize, each=2)
y_gill$samples
```

    ##              group lib.size norm.factors
    ## 14B-G_S25-Me     1   916733            1
    ## 14B-G_S25-Un     1   916733            1
    ## 15G-G_S24-Me     1   306905            1
    ## 15G-G_S24-Un     1   306905            1
    ## 18W-G_S41-Me     1   286880            1
    ## 18W-G_S41-Un     1   286880            1
    ## 21B-G_S7-Me      1  1762059            1
    ## 21B-G_S7-Un      1  1762059            1
    ## 29W-G_S76-Me     1   586501            1
    ## 29W-G_S76-Un     1   586501            1
    ## 29Y-G_S68-Me     1   314022            1
    ## 29Y-G_S68-Un     1   314022            1
    ## 2B-G_S17-Me      1  1286219            1
    ## 2B-G_S17-Un      1  1286219            1
    ## 40Y-G_S63-Me     1   871450            1
    ## 40Y-G_S63-Un     1   871450            1
    ## 45G-G_S69-Me     1   444147            1
    ## 45G-G_S69-Un     1   444147            1
    ## 47W-G_S80-Me     1   798462            1
    ## 47W-G_S80-Un     1   798462            1
    ## 50B-G_S83-Me     1  1312082            1
    ## 50B-G_S83-Un     1  1312082            1
    ## 51B-G_S18-Me     1  1277367            1
    ## 51B-G_S18-Un     1  1277367            1
    ## 55B-G_S40-Me     1   782845            1
    ## 55B-G_S40-Un     1   782845            1
    ## 55Y-G_S10-Me     1   374541            1
    ## 55Y-G_S10-Un     1   374541            1
    ## 56Y-G_S53-Me     1   681562            1
    ## 56Y-G_S53-Un     1   681562            1
    ## 67B-G_S8-Me      1  1889474            1
    ## 67B-G_S8-Un      1  1889474            1
    ## 67W-G_S86-Me     1   975054            1
    ## 67W-G_S86-Un     1   975054            1
    ## 75Y-G_S9-Me      1   595830            1
    ## 75Y-G_S9-Un      1   595830            1
    ## 76Y-G_S85-Me     1   670282            1
    ## 76Y-G_S85-Un     1   670282            1
    ## 90W-G_S67-Me     1   361040            1
    ## 90W-G_S67-Un     1   361040            1

``` r
### Compute the corresponding methylation summary from the methylated and unmethylated counts ###
Me_gill <- y_gill$counts[, Methylation=="Me"]
Un_gill <- y_gill$counts[, Methylation=="Un"]

### Calculating postfiltered genome-wide coverage 
# Coverage_postfilter<-Un_gill + Me_gill
# Coverage_postfilter[Coverage_postfilter == 0] <- NA
# mean_coverage_postfilter <- colMeans(Coverage_postfilter, na.rm = TRUE)
# mean_coverage_postfilter 

### Calculating a methylation proportion matrix
prop_meth_matrix_gill <- Me_gill/(Me_gill+Un_gill)
#Postfiltered total # of CpGs: colSums(!is.na(prop_meth_matrix_gill))

# ### Calculating average genome wide methylation percentage postfiltering
# colMeans(prop_meth_matrix_gill,na.rm=TRUE)*100
```

## Define original and alternative gill models

``` r
### We want to use the exposed site as reference level for both transplant and origin site effect analyses, need to relevel for transplant site effect ###
meta_data_gill$transplant_site <- relevel(factor(meta_data_gill$transplant_site), ref = "exposed")
meta_data_gill$origin_site <- relevel(factor(meta_data_gill$origin_site), ref = "exposed")

### Create a design matrix, using origin site, transplant site, and final shell length as fixed effects ###
designSL_gill <- model.matrix(~0 + origin_site + transplant_site +
                           LENGTH_INIT..mm., 
                           data=meta_data_gill)

# Initial


# Initial plus growth rate

# Merge with growth rate metadata
growth_data_gill <- merge(meta_data_gill,
                          data.frame(
                            meth = shell_df$sample.ID,
                            g_rate = shell_df$mm_per_d,
                            orig_init = shell_df$LENGTH_INIT..mm.,
                            orig_final = shell_df$LENGTH_FINAL..mm.),
                          by = "meth")

designSL_gill_v2 <- model.matrix(~0 + origin_site + transplant_site +
                           LENGTH_INIT..mm. + g_rate, 
                           data = growth_data_gill)

# Just growth rate
designSL_gill_v3 <- model.matrix(~0 + origin_site + transplant_site +
                           g_rate, 
                           data = growth_data_gill)

### Expand to the full design matrix modeling the sample and methylation effects ###
design_gill <- modelMatrixMeth(designSL_gill)
design_gillv2 <- modelMatrixMeth(designSL_gill_v2)
design_gillv3 <- modelMatrixMeth(designSL_gill_v3)

### Dispersion estimation ###
y_gill <- estimateDisp(y_gill, design = design_gill, robust = TRUE)
y_gill2 <- estimateDisp(y_gill, design = design_gillv2, robust = TRUE)
y_gil3 <- estimateDisp(y_gill, design = design_gillv3, robust = TRUE)

### Create the BCV plot for visualization ### 
plotBCV(y_gill)
```

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotBCV(y_gill2)
```

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
plotBCV(y_gil3)
```

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
### Testing for differentially methylated CpG loci ###

### fit NB GLMs for all the CpG loci using the glmFit function in edgeR ###
fit_gill <- glmFit(y_gill, design_gill)
fit_gill2 <- glmFit(y_gill2, design_gillv2)
fit_gill3 <- glmFit(y_gil3, design_gillv3)

### Testing for differentially methylated CpG sites between different treatment groups using likelihood ratio tests ###

### Origin site effects ###
contr_origin_gill <- makeContrasts(Origin = origin_siteprotected-origin_siteexposed, levels = design_gill)
contr_origin_gill2 <- makeContrasts(Origin = origin_siteprotected-origin_siteexposed, levels = design_gillv2)
contr_origin_gill3 <- makeContrasts(Origin = origin_siteprotected-origin_siteexposed, levels = design_gillv3)

lrt_origin_gill <- glmLRT(fit_gill, contrast=contr_origin_gill)
lrt_origin_gill2 <- glmLRT(fit_gill2, contrast=contr_origin_gill2)
lrt_origin_gill3 <- glmLRT(fit_gill3, contrast=contr_origin_gill3)

summary(decideTests(lrt_origin_gill))
```

    ##        -1*origin_siteexposed 1*origin_siteprotected
    ## Down                                             55
    ## NotSig                                        76528
    ## Up                                               53

``` r
summary(decideTests(lrt_origin_gill2))
```

    ##        -1*origin_siteexposed 1*origin_siteprotected
    ## Down                                             64
    ## NotSig                                        76523
    ## Up                                               49

``` r
summary(decideTests(lrt_origin_gill3))
```

    ##        -1*origin_siteexposed 1*origin_siteprotected
    ## Down                                             43
    ## NotSig                                        76545
    ## Up                                               48

``` r
### Transplant site effects (exposed is the reference level) ###
contr_trans_gill <- makeContrasts(Transplant = transplant_siteprotected, levels = design_gill)
contr_trans_gill2 <- makeContrasts(Transplant = transplant_siteprotected, levels = design_gillv2)
contr_trans_gill3 <- makeContrasts(Transplant = transplant_siteprotected, levels = design_gillv3)

lrt_trans_gill <- glmLRT(fit_gill, contrast=contr_trans_gill)
lrt_trans_gill2 <- glmLRT(fit_gill2, contrast=contr_trans_gill2)
lrt_trans_gill3 <- glmLRT(fit_gill3, contrast=contr_trans_gill3)

summary(decideTests(lrt_trans_gill))
```

    ##        1*transplant_siteprotected
    ## Down                           67
    ## NotSig                      76400
    ## Up                            169

``` r
summary(decideTests(lrt_trans_gill2))
```

    ##        1*transplant_siteprotected
    ## Down                           28
    ## NotSig                      76553
    ## Up                             55

``` r
summary(decideTests(lrt_trans_gill3))
```

    ##        1*transplant_siteprotected
    ## Down                           25
    ## NotSig                      76562
    ## Up                             49

``` r
### Length effect ###
contr_length_gill <- makeContrasts(Length = LENGTH_INIT..mm., levels = design_gill)

contr_length_gill2 <- makeContrasts(Length = LENGTH_INIT..mm., levels = design_gillv2)
contr_grate_gill2 <- makeContrasts(grate = g_rate, levels = design_gillv2)

contr_grate_gill3 <- makeContrasts(grate = g_rate, levels = design_gillv3)

lrt_length_gill <- glmLRT(fit_gill, contrast=contr_length_gill)

lrt_length_gill2 <- glmLRT(fit_gill2, contrast=contr_length_gill2)
lrt_grate_gill2 <- glmLRT(fit_gill2, contrast=contr_grate_gill2)

lrt_grate_gill3 <- glmLRT(fit_gill3, contrast=contr_grate_gill3)

### Wrangle data for generating volano plot of CpG DM ###
### Correct p-values using BH method ###
lrt_origin_gill$table$FDR <- p.adjust(lrt_origin_gill$table$PValue, method = "BH")
lrt_origin_gill2$table$FDR <- p.adjust(lrt_origin_gill2$table$PValue, method = "BH")
lrt_origin_gill3$table$FDR <- p.adjust(lrt_origin_gill3$table$PValue, method = "BH")

lrt_trans_gill$table$FDR <- p.adjust(lrt_trans_gill$table$PValue, method = "BH")
lrt_trans_gill2$table$FDR <- p.adjust(lrt_trans_gill2$table$PValue, method = "BH")
lrt_trans_gill3$table$FDR <- p.adjust(lrt_trans_gill3$table$PValue, method = "BH")

lrt_length_gill$table$FDR <- p.adjust(lrt_length_gill$table$PValue, method = "BH")

lrt_length_gill2$table$FDR <- p.adjust(lrt_length_gill2$table$PValue, method = "BH")
lrt_grate_gill2$table$FDR <- p.adjust(lrt_grate_gill2$table$PValue, method = "BH")

lrt_grate_gill3$table$FDR <- p.adjust(lrt_grate_gill3$table$PValue, method = "BH")
```

## Use DM info to evaluate which predictors explain more DM

``` r
# If growth rate is a useful additional parameter, more DM CpGs should associate with growth in v3 length than len in v1
# If addition of growth reveals length DM, # of len DM CpGs should increase
# Init len only
summary(decideTests(lrt_length_gill))
```

    ##        1*LENGTH_INIT..mm.
    ## Down                   98
    ## NotSig              76478
    ## Up                     60

``` r
98+60
```

    ## [1] 158

``` r
# Full len DM CpGs
summary(decideTests(lrt_length_gill2))
```

    ##        1*LENGTH_INIT..mm.
    ## Down                   38
    ## NotSig              76537
    ## Up                     61

``` r
38+61 # 33% reduction in DM
```

    ## [1] 99

``` r
# V2
summary(decideTests(lrt_grate_gill2))
```

    ##        1*g_rate
    ## Down         50
    ## NotSig    76523
    ## Up           63

``` r
50+63+38+61 # Expected increase in total DM compared to length alone
```

    ## [1] 212

``` r
summary(decideTests(lrt_grate_gill3))
```

    ##        1*g_rate
    ## Down         79
    ## NotSig    76496
    ## Up           61

``` r
79+61 # 140...  fewer than 158
```

    ## [1] 140

``` r
# How many of these CpGs are similar? Probably very few
table(
  row.names(
  filter(
    as.data.frame(decideTests(lrt_grate_gill3)), `1*g_rate` %in% c(1,-1)
    )
  ) %in% row.names(
  filter(
    as.data.frame(decideTests(lrt_length_gill)), `1*LENGTH_INIT..mm.` %in% c(1,-1)
    )
  )
)
```

    ## 
    ## FALSE  TRUE 
    ##    99    41

``` r
# How much overlap? 31%... fairly high
41/(91+41)
```

    ## [1] 0.3106061

Now, let’s calc AIC for a few of these top transplant genes

``` r
library(MASS)

## Initial length vs length + growth

# Identify 158 significant DM CpGs from gill transplant effect (model v1)
top158_ids <- rownames(
  lrt_trans_gill$table[order(lrt_trans_gill$table$PValue), ][1:158, ]
)

# Subset to top 10 CpGs
Me_top  <- Me_gill[top158_ids, ]
Un_top  <- Un_gill[top158_ids, ]
Cov_top <- Me_top + Un_top

# Calc model AIC
aic_results <- sapply(1:length(top158_ids), function(i) {
  
  df <- data.frame(
    Me = as.numeric(Me_top[i, ]),
    Cov = as.numeric(Cov_top[i, ]),
    origin_site = meta_data_gill$origin_site,
    transplant_site = meta_data_gill$transplant_site,
    LENGTH_INIT = meta_data_gill$LENGTH_INIT..mm.
  )
  
  # Remove 0 cov rows
  df <- df[df$Cov > 0, ]
  
  fit <- glm.nb(
    Me ~ 0 + origin_site + transplant_site + LENGTH_INIT + 
      offset(log(Cov)),
    data = df
  )
  
  AIC(fit)
})

names(aic_results) <- top158_ids
aic_results <- as.data.frame(aic_results)

# V2 AIC results
aic_results2 <- sapply(1:length(top158_ids), function(i) {
  
  df <- data.frame(
    Me = as.numeric(Me_top[i, ]),
    Cov = as.numeric(Cov_top[i, ]),
    origin_site = meta_data_gill$origin_site,
    transplant_site = meta_data_gill$transplant_site,
    LENGTH_INIT = meta_data_gill$LENGTH_INIT..mm.,
    g_rate = growth_data_gill$g_rate
  )
  
  # Remove zero coverage rows
  df <- df[df$Cov > 0, ]
  
  fit <- glm.nb(
    Me ~ 0 + origin_site + transplant_site + LENGTH_INIT + g_rate +
      offset(log(Cov)),
    data = df
  )
  
  AIC(fit)
})

aic_results2 <- as.data.frame(aic_results2)

aic_comp <- data.frame(CpG = row.names(aic_results),
                       AIC_len = aic_results$aic_results,
                       AIC_len_grate = aic_results2$aic_results2)

aic_comp$aic_simp_vs_full <- aic_comp$AIC_len - aic_comp$AIC_len_grate

# Plot AIC difference
ggplot(data = aic_comp,
       aes(x = aic_simp_vs_full)) +
  geom_histogram(binwidth = .5, fill = "grey", color = "black") +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic() +
  labs(x = "AIC(len) - AIC(len + growth)", y = "Count (CpG models)")
```

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Identify 10 most-significant DM CpGs from gill transplant v1
top158_ids <- rownames(
  lrt_trans_gill$table[order(lrt_trans_gill$table$PValue), ][1:158, ]
)



## Initial growth rate vs. growth

# V3 AIC results
aic_results3 <- sapply(1:length(top158_ids), function(i) {
  
  df <- data.frame(
    Me = as.numeric(Me_top[i, ]),
    Cov = as.numeric(Cov_top[i, ]),
    origin_site = meta_data_gill$origin_site,
    transplant_site = meta_data_gill$transplant_site,
    LENGTH_INIT = meta_data_gill$LENGTH_INIT..mm.,
    g_rate = growth_data_gill$g_rate
  )
  
  # Remove zero coverage rows
  df <- df[df$Cov > 0, ]
  
  fit <- glm.nb(
    Me ~ 0 + origin_site + transplant_site + g_rate +
      offset(log(Cov)),
    data = df
  )
  
  AIC(fit)
})

aic_results3 <- as.data.frame(aic_results3)

aic_comp2 <- data.frame(CpG = row.names(aic_results),
                       AIC_len = aic_results$aic_results,
                       AIC_grate = aic_results3$aic_results3)

aic_comp2$aic_simp_vs_full <- aic_comp2$AIC_len - aic_comp2$AIC_grate

# Plot AIC difference
ggplot(data = aic_comp2,
       aes(x = aic_simp_vs_full)) +
  geom_histogram(binwidth = .5, fill = "grey", color = "black") +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic() +
  labs(x = "AIC(len) - AIC(growth)", y = "Count (CpG models)")
```

![](sb_growth_reanalysis_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->
