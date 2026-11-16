# aRMS_SS_DESeq2


## aRMS v SS

This notebook will be analyzing the expression profiles of synovial
sarcoma (SS) and alveolar rhabdomyosarcoma (aRMS) tumors from the v25
riboD and polyA compendia.

The goal is to determine the extent that genes get called as
differentially expressed as a result of the transcript enrichment method
of the library preparation process.

The comparisons we’ll be making (SS = synovial sarcoma, aRMS = alveolar
rhabdomyosarcoma):

- aRMS polyA vs SS polyA

- aRMS riboD vs SS riboD

- aRMS polyA vs SS riboD

- aRMS riboD vs SS polyA

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ✔ forcats   1.0.1     ✔ stringr   1.5.2
    ✔ ggplot2   4.0.0     ✔ tibble    3.3.0
    ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ✔ purrr     1.1.0     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(DESeq2)
```

    Loading required package: S4Vectors
    Loading required package: stats4
    Loading required package: BiocGenerics
    Loading required package: generics

    Attaching package: 'generics'

    The following object is masked from 'package:lubridate':

        as.difftime

    The following object is masked from 'package:dplyr':

        explain

    The following objects are masked from 'package:base':

        as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
        setequal, union


    Attaching package: 'BiocGenerics'

    The following object is masked from 'package:dplyr':

        combine

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
        unsplit, which.max, which.min


    Attaching package: 'S4Vectors'

    The following objects are masked from 'package:lubridate':

        second, second<-

    The following objects are masked from 'package:dplyr':

        first, rename

    The following object is masked from 'package:tidyr':

        expand

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Attaching package: 'IRanges'

    The following object is masked from 'package:lubridate':

        %within%

    The following objects are masked from 'package:dplyr':

        collapse, desc, slice

    The following object is masked from 'package:purrr':

        reduce

    Loading required package: GenomicRanges
    Loading required package: Seqinfo
    Loading required package: SummarizedExperiment
    Loading required package: MatrixGenerics
    Loading required package: matrixStats

    Attaching package: 'matrixStats'

    The following object is masked from 'package:dplyr':

        count


    Attaching package: 'MatrixGenerics'

    The following objects are masked from 'package:matrixStats':

        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars

    Loading required package: Biobase
    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.


    Attaching package: 'Biobase'

    The following object is masked from 'package:MatrixGenerics':

        rowMedians

    The following objects are masked from 'package:matrixStats':

        anyMissing, rowMedians

``` r
# sample metadata
# This includes the top 15 samples from SS polyA, SS riboD, aRMS polyA, aRMS riboD
# The top 15 samples came from the correlation analysis, where this tsv file was generated
filtered_SS_aRMS_list <- read_tsv("../../input_data/filtered_SS_aRMS_list.tsv", show_col_types = FALSE)

expression_SS_aRMS <- read_tsv("../../input_data/expression_SS_aRMS.tsv.gz", show_col_types = FALSE)

gene_names <- read_tsv("../../input_data/EnsGeneID_Hugo_Observed_Conversions.txt", show_col_types = FALSE)
```

That expression file contains all expression information for all SS and
aRMS samples across the compendia (170 samples). We want only the
expression information for the samples we’re interested in, based on the
correlation analysis we did previously (60 samples).

``` r
# selecting only the 15 most correlated samples for each category
interested_expression <- expression_SS_aRMS %>%
  select(all_of(c("Gene", filtered_SS_aRMS_list$term)))
```

Since filtering out genes with a count less than 10 across all samples
will filter out relevant genes, we’ll filter out genes with a count less
than 10 in 11 or more samples (i.e., we’re keeping genes that have a
count greater than 10 in at least 10 samples).

For example, if Gene A is expressed in SS and not in aRMS, but it’s only
picked up by polyA-selection, then only SS-polyA would have expression
for it, and the gene would be filtered out in the pre-group filtering
step below. However, we would want to capture that gene information,
because it is both biologically differentially expressed (high
expression in SS, low expression in aRMS) and represents an artifact of
library prep (high expression in SS-polyA, low expression in SS-riboD).

Thus, we filter for genes with expression above 10 in at least 10
samples, which would capture the expression of Gene A in the example
above.

``` r
# Explanation:
# across(where(is.numeric)) is selecting all the numeric columns
# i.e., all the sample columns. This deselects the Gene column, since that is a character column

# rowSums() counts how many samples meet the threshold for each gene,
# i.e., have a value of at least 10

# then keeps only the genes that meet that threshold in at least 10 samples
updated_expression <- interested_expression %>%
  filter(rowSums(across(where(is.numeric)) >= 10) >= 10)
```

### aRMS polyA vs SS polyA

``` r
# metadata
aRMS_polyA_SS_polyA_meta <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep == "aRMS_polyA" | disease_and_prep == "SS_polyA")
# the pipe | means we're taking all samples called SS_riboD OR SS_polyA in column disease_and_prep

# expression
aRMS_polyA_SS_polyA_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", aRMS_polyA_SS_polyA_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_aRMS_polyA_SS_polyA_updated_expr <- aRMS_polyA_SS_polyA_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_aRMS_polyA_SS_polyA_meta <- aRMS_polyA_SS_polyA_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_aRMS_polyA_SS_polyA_meta) %in% colnames(fix_aRMS_polyA_SS_polyA_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_aRMS_polyA_SS_polyA_meta) == colnames(fix_aRMS_polyA_SS_polyA_updated_expr))
```

    [1] TRUE

``` r
# format data
round_aRMS_polyA_SS_polyA_updated_expr = round(fix_aRMS_polyA_SS_polyA_updated_expr)

aRMS_polyA_SS_polyA_updated_expr_matrix <- as.matrix(round_aRMS_polyA_SS_polyA_updated_expr)

fix_aRMS_polyA_SS_polyA_meta$disease_and_prep <- as.factor(fix_aRMS_polyA_SS_polyA_meta$disease_and_prep)
```

``` r
# set reference levels
# we want to make sure the references are consistent across the comparisons
levels(fix_aRMS_polyA_SS_polyA_meta$disease_and_prep)
```

    [1] "aRMS_polyA" "SS_polyA"  

``` r
# set up DESeq function
dds_aRMS_polyA_SS_polyA_updated <- DESeqDataSetFromMatrix(countData = aRMS_polyA_SS_polyA_updated_expr_matrix,
                              colData = fix_aRMS_polyA_SS_polyA_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
dds_aRMS_polyA_SS_polyA_updated
```

    class: DESeqDataSet 
    dim: 27327 30 
    metadata(1): version
    assays(1): counts
    rownames(27327): ENSG00000000003.14 ENSG00000000005.5 ...
      ENSG00000282787.1 ENSG00000282798.1
    rowData names(0):
    colnames(30): THR39_1373_S01 TCGA-WK-A8XT-01 ... THR29_0738_S01
      THR29_0720_S01
    colData names(1): disease_and_prep

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
aRMS_polyA_SS_polyA_updated_genes_to_keep <- rowSums(counts(dds_aRMS_polyA_SS_polyA_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_aRMS_polyA_SS_polyA_updated <- dds_aRMS_polyA_SS_polyA_updated[which(aRMS_polyA_SS_polyA_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_aRMS_polyA_SS_polyA_updated_object <- DESeq(dds_aRMS_polyA_SS_polyA_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 749 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
result_dds_aRMS_polyA_SS_polyA_updated <- results(dds_aRMS_polyA_SS_polyA_updated_object)
# summarize results table
summary(result_dds_aRMS_polyA_SS_polyA_updated, alpha = 0.05)
```


    out of 26946 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 7659, 28%
    LFC < 0 (down)     : 9653, 36%
    outliers [1]       : 0, 0%
    low counts [2]     : 3, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# convert results to dataframe
result_dds_aRMS_polyA_SS_polyA_updated_df <- result_dds_aRMS_polyA_SS_polyA_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
aRMS_polyA_SS_polyA_LFC1great_padj <- result_dds_aRMS_polyA_SS_polyA_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
aRMS_polyA_SS_polyA_LFC1great_padj
```

    [1] 4622

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
aRMS_polyA_SS_polyA_LFC1less_padj <- result_dds_aRMS_polyA_SS_polyA_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
aRMS_polyA_SS_polyA_LFC1less_padj
```

    [1] 6759

``` r
# convert EnsemblID to HugoID for gene name readability
hugo_results_aRMS_polyA_SS_polyA_updated <- result_dds_aRMS_polyA_SS_polyA_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_aRMS_polyA_SS_polyA_updated, path = "../../output_data/aRMS_SS/hugo_results_aRMS_polyA_SS_polyA_updated.tsv.gz")
```

    Warning: The `path` argument of `write_tsv()` is deprecated as of readr 1.4.0.
    ℹ Please use the `file` argument instead.

### aRMS riboD vs SS riboD

``` r
# metadata
aRMS_riboD_SS_riboD_meta <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep == "aRMS_riboD" | disease_and_prep == "SS_riboD")

# expression
aRMS_riboD_SS_riboD_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", aRMS_riboD_SS_riboD_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_aRMS_riboD_SS_riboD_updated_expr <- aRMS_riboD_SS_riboD_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_aRMS_riboD_SS_riboD_meta <- aRMS_riboD_SS_riboD_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_aRMS_riboD_SS_riboD_meta) %in% colnames(fix_aRMS_riboD_SS_riboD_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_aRMS_riboD_SS_riboD_meta) == colnames(fix_aRMS_riboD_SS_riboD_updated_expr))
```

    [1] TRUE

``` r
# format data
round_aRMS_riboD_SS_riboD_updated_expr = round(fix_aRMS_riboD_SS_riboD_updated_expr)

aRMS_riboD_SS_riboD_updated_expr_matrix <- as.matrix(round_aRMS_riboD_SS_riboD_updated_expr)

fix_aRMS_riboD_SS_riboD_meta$disease_and_prep <- as.factor(fix_aRMS_riboD_SS_riboD_meta$disease_and_prep)
```

``` r
# set reference levels
# we want to make sure the references are consistent across the comparisons
levels(fix_aRMS_riboD_SS_riboD_meta$disease_and_prep)
```

    [1] "aRMS_riboD" "SS_riboD"  

``` r
# set up DESeq function
dds_aRMS_riboD_SS_riboD_updated <- DESeqDataSetFromMatrix(countData = aRMS_riboD_SS_riboD_updated_expr_matrix,
                              colData = fix_aRMS_riboD_SS_riboD_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
aRMS_riboD_SS_riboD_updated_genes_to_keep <- rowSums(counts(dds_aRMS_riboD_SS_riboD_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_aRMS_riboD_SS_riboD_updated <- dds_aRMS_riboD_SS_riboD_updated[which(aRMS_riboD_SS_riboD_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_aRMS_riboD_SS_riboD_updated_object <- DESeq(dds_aRMS_riboD_SS_riboD_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 1651 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
res_aRMS_riboD_SS_riboD_updated <- results(dds_aRMS_riboD_SS_riboD_updated_object)

# summarize results table
summary(res_aRMS_riboD_SS_riboD_updated, alpha = 0.05)
```


    out of 27309 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 6541, 24%
    LFC < 0 (down)     : 6594, 24%
    outliers [1]       : 0, 0%
    low counts [2]     : 3, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# first we convert to df
res_aRMS_riboD_SS_riboD_updated_df <- res_aRMS_riboD_SS_riboD_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
aRMS_riboD_SS_riboD_LFC1great_padj <- res_aRMS_riboD_SS_riboD_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
aRMS_riboD_SS_riboD_LFC1great_padj
```

    [1] 3921

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
aRMS_riboD_SS_riboD_LFC1less_padj <- res_aRMS_riboD_SS_riboD_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
aRMS_riboD_SS_riboD_LFC1less_padj
```

    [1] 4981

``` r
hugo_results_aRMS_riboD_SS_riboD_updated <- res_aRMS_riboD_SS_riboD_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_aRMS_riboD_SS_riboD_updated, path = "../../output_data/aRMS_SS/hugo_results_aRMS_riboD_SS_riboD_updated.tsv.gz")
```

### **aRMS polyA vs SS riboD**

``` r
# metadata
aRMS_polyA_SS_riboD_meta <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep == "aRMS_polyA" | disease_and_prep == "SS_riboD")
# the pipe | means we're taking all samples called SS_riboD OR SS_polyA in column disease_and_prep

# expression
aRMS_polyA_SS_riboD_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", aRMS_polyA_SS_riboD_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_aRMS_polyA_SS_riboD_updated_expr <- aRMS_polyA_SS_riboD_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_aRMS_polyA_SS_riboD_meta <- aRMS_polyA_SS_riboD_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_aRMS_polyA_SS_riboD_meta) %in% colnames(fix_aRMS_polyA_SS_riboD_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_aRMS_polyA_SS_riboD_meta) == colnames(fix_aRMS_polyA_SS_riboD_updated_expr))
```

    [1] TRUE

``` r
# format data
round_aRMS_polyA_SS_riboD_updated_expr = round(fix_aRMS_polyA_SS_riboD_updated_expr)

aRMS_polyA_SS_riboD_updated_expr_matrix <- as.matrix(round_aRMS_polyA_SS_riboD_updated_expr)

fix_aRMS_polyA_SS_riboD_meta$disease_and_prep <- as.factor(fix_aRMS_polyA_SS_riboD_meta$disease_and_prep)
```

``` r
# check reference levels
levels(fix_aRMS_polyA_SS_riboD_meta$disease_and_prep)
```

    [1] "aRMS_polyA" "SS_riboD"  

``` r
# set up DESeq function
dds_aRMS_polyA_SS_riboD_updated <- DESeqDataSetFromMatrix(countData = aRMS_polyA_SS_riboD_updated_expr_matrix,
                              colData = fix_aRMS_polyA_SS_riboD_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
dds_aRMS_polyA_SS_riboD_updated
```

    class: DESeqDataSet 
    dim: 27327 30 
    metadata(1): version
    assays(1): counts
    rownames(27327): ENSG00000000003.14 ENSG00000000005.5 ...
      ENSG00000282787.1 ENSG00000282798.1
    rowData names(0):
    colnames(30): THR51_4556_S01 THR51_4558_S01 ... THR29_0738_S01
      THR29_0720_S01
    colData names(1): disease_and_prep

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
aRMS_polyA_SS_riboD_updated_genes_to_keep <- rowSums(counts(dds_aRMS_polyA_SS_riboD_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_aRMS_polyA_SS_riboD_updated <- dds_aRMS_polyA_SS_riboD_updated[which(aRMS_polyA_SS_riboD_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_aRMS_polyA_SS_riboD_updated_object <- DESeq(dds_aRMS_polyA_SS_riboD_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 1277 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
result_dds_aRMS_polyA_SS_riboD_updated <- results(dds_aRMS_polyA_SS_riboD_updated_object)
# summarize results table
summary(result_dds_aRMS_polyA_SS_riboD_updated, alpha = 0.05)
```


    out of 27325 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 11310, 41%
    LFC < 0 (down)     : 8312, 30%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# convert results to dataframe
result_dds_aRMS_polyA_SS_riboD_updated_df <- result_dds_aRMS_polyA_SS_riboD_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
aRMS_polyA_SS_riboD_LFC1great_padj <- result_dds_aRMS_polyA_SS_riboD_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
aRMS_polyA_SS_riboD_LFC1great_padj
```

    [1] 8951

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
aRMS_polyA_SS_riboD_LFC1less_padj <- result_dds_aRMS_polyA_SS_riboD_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
aRMS_polyA_SS_riboD_LFC1less_padj
```

    [1] 5711

``` r
hugo_results_aRMS_polyA_SS_riboD_updated <- result_dds_aRMS_polyA_SS_riboD_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_aRMS_polyA_SS_riboD_updated, path = "../../output_data/aRMS_SS/hugo_results_aRMS_polyA_SS_riboD_updated.tsv.gz")
```

### **aRMS riboD vs SS polyA**

``` r
# metadata
aRMS_riboD_SS_polyA_meta <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep == "aRMS_riboD" | disease_and_prep == "SS_polyA")

# expression
aRMS_riboD_SS_polyA_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", aRMS_riboD_SS_polyA_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_aRMS_riboD_SS_polyA_updated_expr <- aRMS_riboD_SS_polyA_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_aRMS_riboD_SS_polyA_meta <- aRMS_riboD_SS_polyA_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_aRMS_riboD_SS_polyA_meta) %in% colnames(fix_aRMS_riboD_SS_polyA_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_aRMS_riboD_SS_polyA_meta) == colnames(fix_aRMS_riboD_SS_polyA_updated_expr))
```

    [1] TRUE

``` r
# format data
round_aRMS_riboD_SS_polyA_updated_expr = round(fix_aRMS_riboD_SS_polyA_updated_expr)

aRMS_riboD_SS_polyA_updated_expr_matrix <- as.matrix(round_aRMS_riboD_SS_polyA_updated_expr)

fix_aRMS_riboD_SS_polyA_meta$disease_and_prep <- as.factor(fix_aRMS_riboD_SS_polyA_meta$disease_and_prep)
```

``` r
# check reference levels
levels(fix_aRMS_riboD_SS_polyA_meta$disease_and_prep)
```

    [1] "aRMS_riboD" "SS_polyA"  

``` r
# set up DESeq function
dds_aRMS_riboD_SS_polyA_updated <- DESeqDataSetFromMatrix(countData = aRMS_riboD_SS_polyA_updated_expr_matrix,
                              colData = fix_aRMS_riboD_SS_polyA_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
aRMS_riboD_SS_polyA_updated_genes_to_keep <- rowSums(counts(dds_aRMS_riboD_SS_polyA_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_aRMS_riboD_SS_polyA_updated <- dds_aRMS_riboD_SS_polyA_updated[which(aRMS_riboD_SS_polyA_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_aRMS_riboD_SS_polyA_updated_object <- DESeq(dds_aRMS_riboD_SS_polyA_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 959 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
res_aRMS_riboD_SS_polyA_updated <- results(dds_aRMS_riboD_SS_polyA_updated_object)

# summarize results table
summary(res_aRMS_riboD_SS_polyA_updated, alpha = 0.05)
```


    out of 27324 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 7870, 29%
    LFC < 0 (down)     : 12556, 46%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# first we convert to df
res_aRMS_riboD_SS_polyA_updated_df <- res_aRMS_riboD_SS_polyA_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
aRMS_riboD_SS_polyA_LFC1great_padj <- res_aRMS_riboD_SS_polyA_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
aRMS_riboD_SS_polyA_LFC1great_padj
```

    [1] 4894

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
aRMS_riboD_SS_polyA_LFC1less_padj <- res_aRMS_riboD_SS_polyA_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
aRMS_riboD_SS_polyA_LFC1less_padj
```

    [1] 10438

``` r
hugo_results_aRMS_riboD_SS_polyA_updated <- res_aRMS_riboD_SS_polyA_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_aRMS_riboD_SS_polyA_updated, path = "../../output_data/aRMS_SS/hugo_results_aRMS_riboD_SS_polyA_updated.tsv.gz")
```

### Session Info

``` r
sessioninfo::session_info()
```

    ─ Session info ───────────────────────────────────────────────────────────────
     setting  value
     version  R version 4.5.2 (2025-10-31)
     os       macOS Tahoe 26.5
     system   aarch64, darwin20
     ui       X11
     language (EN)
     collate  en_US.UTF-8
     ctype    en_US.UTF-8
     tz       America/Los_Angeles
     date     2026-06-29
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     ! package              * version  date (UTC) lib source
     P abind                  1.4-8    2024-09-12 [?] CRAN (R 4.5.0)
     P Biobase              * 2.70.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P BiocGenerics         * 0.56.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P BiocManager            1.30.26  2025-06-05 [?] CRAN (R 4.5.0)
     P BiocParallel           1.44.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P bit                    4.6.0    2025-03-06 [?] RSPM
     P bit64                  4.6.0-1  2025-01-16 [?] CRAN (R 4.5.0)
     P cli                    3.6.5    2025-04-23 [?] CRAN (R 4.5.0)
     P codetools              0.2-20   2024-03-31 [?] CRAN (R 4.5.2)
     P crayon                 1.5.3    2024-06-20 [?] RSPM
     P DelayedArray           0.36.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P DESeq2               * 1.50.2   2025-11-12 [?] https://bioc-release.r-universe.dev (R 4.5.2)
     P digest                 0.6.37   2024-08-19 [?] CRAN (R 4.5.0)
     P dplyr                * 1.1.4    2023-11-17 [?] CRAN (R 4.5.0)
     P evaluate               1.0.5    2025-08-27 [?] RSPM
     P farver                 2.1.2    2024-05-13 [?] RSPM
     P fastmap                1.2.0    2024-05-15 [?] RSPM
     P forcats              * 1.0.1    2025-09-25 [?] RSPM
     P generics             * 0.1.4    2025-05-09 [?] RSPM
     P GenomicRanges        * 1.62.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P ggplot2              * 4.0.0    2025-09-11 [?] CRAN (R 4.5.0)
     P glue                   1.8.0    2024-09-30 [?] CRAN (R 4.5.0)
     P gtable                 0.3.6    2024-10-25 [?] RSPM
     P hms                    1.1.4    2025-10-17 [?] RSPM
     P htmltools              0.5.8.1  2024-04-04 [?] CRAN (R 4.5.0)
     P IRanges              * 2.44.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P jsonlite               2.0.0    2025-03-27 [?] RSPM
     P knitr                  1.50     2025-03-16 [?] CRAN (R 4.5.0)
     P lattice                0.22-7   2025-04-02 [?] CRAN (R 4.5.2)
     P lifecycle              1.0.4    2023-11-07 [?] CRAN (R 4.5.0)
     P locfit                 1.5-9.12 2025-03-05 [?] CRAN (R 4.5.0)
     P lubridate            * 1.9.4    2024-12-08 [?] CRAN (R 4.5.0)
     P magrittr               2.0.4    2025-09-12 [?] CRAN (R 4.5.0)
     P Matrix                 1.7-4    2025-08-28 [?] CRAN (R 4.5.2)
     P MatrixGenerics       * 1.22.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P matrixStats          * 1.5.0    2025-01-07 [?] CRAN (R 4.5.0)
     P pillar                 1.11.1   2025-09-17 [?] RSPM
     P pkgconfig              2.0.3    2019-09-22 [?] RSPM
     P purrr                * 1.1.0    2025-07-10 [?] CRAN (R 4.5.0)
     P R6                     2.6.1    2025-02-15 [?] RSPM
     P RColorBrewer           1.1-3    2022-04-03 [?] RSPM
     P Rcpp                   1.1.0    2025-07-02 [?] CRAN (R 4.5.0)
     P readr                * 2.1.5    2024-01-10 [?] CRAN (R 4.5.0)
       renv                   1.1.5    2025-07-24 [1] CRAN (R 4.5.0)
     P rlang                  1.2.0    2026-04-06 [?] RSPM
     P rmarkdown              2.30     2025-09-28 [?] CRAN (R 4.5.0)
     P rstudioapi             0.17.1   2024-10-22 [?] CRAN (R 4.5.0)
     P S4Arrays               1.10.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P S4Vectors            * 0.48.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P S7                     0.2.0    2024-11-07 [?] CRAN (R 4.5.0)
     P scales                 1.4.0    2025-04-24 [?] RSPM
     P Seqinfo              * 1.0.0    2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P sessioninfo            1.2.3    2025-02-05 [?] CRAN (R 4.5.0)
     P SparseArray            1.10.1   2025-10-31 [?] Bioconductor 3.22 (R 4.5.1)
     P stringi                1.8.7    2025-03-27 [?] RSPM
     P stringr              * 1.5.2    2025-09-08 [?] CRAN (R 4.5.0)
     P SummarizedExperiment * 1.40.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P tibble               * 3.3.0    2025-06-08 [?] CRAN (R 4.5.0)
     P tidyr                * 1.3.1    2024-01-24 [?] CRAN (R 4.5.0)
     P tidyselect             1.2.1    2024-03-11 [?] RSPM
     P tidyverse            * 2.0.0    2023-02-22 [?] RSPM
     P timechange             0.3.0    2024-01-18 [?] CRAN (R 4.5.0)
     P tzdb                   0.5.0    2025-03-15 [?] RSPM
     P vctrs                  0.6.5    2023-12-01 [?] CRAN (R 4.5.0)
     P vroom                  1.6.6    2025-09-19 [?] CRAN (R 4.5.0)
     P withr                  3.0.2    2024-10-28 [?] CRAN (R 4.5.0)
     P xfun                   0.55     2025-12-16 [?] CRAN (R 4.5.2)
     P XVector                0.50.0   2025-10-29 [?] Bioconductor 3.22 (R 4.5.1)
     P yaml                   2.3.10   2024-07-26 [?] CRAN (R 4.5.0)

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/Fig_2/aRMS_SS/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Users/maryke/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
