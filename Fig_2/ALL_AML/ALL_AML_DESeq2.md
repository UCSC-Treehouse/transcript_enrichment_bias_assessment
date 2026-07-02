# ALL_AML_DESeq2


## ALL v AML

This notebook is analyzing the expression profiles of acute myeloid
leukemia (AML) and acute lymphoblastic leukemia (ALL) tumors from the
v25 riboD and polyA compendia.

Goal: determine the extent that genes get called as differentially
expressed as a result of the transcript enrichment method of the library
preparation process.

Use: subsetted list of tumors for each disease in each library prep
method (based on correlation analysis and tumor map analysis done
before)

The comparisons we’ll be making (ALL = acute lymphoblastic leukemia, AML
= acute myeloid leukemia):

- AML polyA vs ALL polyA

- AML riboD vs ALL riboD

- AML polyA vs ALL riboD

- AML riboD vs ALL polyA

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
# This includes the top 20 samples from ALL polyA, ALL riboD, AML polyA, AML riboD
# The top 20 samples came from the correlation analysis, where this tsv file was generated
filtered_ALL_AML_list <- read_tsv("../../input_data/filtered_ALL_AML_list.tsv", show_col_types = FALSE)

expression_ALL_AML <- read_tsv("../../input_data/expression_ALL_AML.tsv.gz", show_col_types = FALSE)

gene_names <- read_tsv("../../input_data/EnsGeneID_Hugo_Observed_Conversions.txt", show_col_types = FALSE)
```

``` r
# selecting only the 20 most correlated samples for each category,
# based on the metadata file
interested_expression <- expression_ALL_AML %>%
  select(all_of(c("Gene", filtered_ALL_AML_list$term)))
```

Since filtering out genes with a count less than 10 across all samples
will filter out relevant genes, we’ll filter out genes with a count less
than 10 in 11 or more samples (i.e., we’re keeping genes that have a
count greater than 10 in at least 10 samples).

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

### **AML polyA vs ALL polyA**

``` r
# metadata
AML_polyA_ALL_polyA_meta <- filtered_ALL_AML_list %>%
  filter(disease_and_prep == "AML_polyA" | disease_and_prep == "ALL_polyA")
# the pipe | means we're taking all samples called ALL_riboD OR ALL_polyA in column disease_and_prep

# expression
AML_polyA_ALL_polyA_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", AML_polyA_ALL_polyA_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_AML_polyA_ALL_polyA_updated_expr <- AML_polyA_ALL_polyA_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_AML_polyA_ALL_polyA_meta <- AML_polyA_ALL_polyA_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_AML_polyA_ALL_polyA_meta) %in% colnames(fix_AML_polyA_ALL_polyA_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_AML_polyA_ALL_polyA_meta) == colnames(fix_AML_polyA_ALL_polyA_updated_expr))
```

    [1] TRUE

``` r
# format data
round_AML_polyA_ALL_polyA_updated_expr = round(fix_AML_polyA_ALL_polyA_updated_expr)

AML_polyA_ALL_polyA_updated_expr_matrix <- as.matrix(round_AML_polyA_ALL_polyA_updated_expr)

fix_AML_polyA_ALL_polyA_meta$disease_and_prep <- as.factor(fix_AML_polyA_ALL_polyA_meta$disease_and_prep)
```

``` r
# set reference levels
# we want to make sure the references are consistent across the comparisons
levels(fix_AML_polyA_ALL_polyA_meta$disease_and_prep)
```

    [1] "ALL_polyA" "AML_polyA"

``` r
# set up DESeq function
dds_AML_polyA_ALL_polyA_updated <- DESeqDataSetFromMatrix(countData = AML_polyA_ALL_polyA_updated_expr_matrix,
                              colData = fix_AML_polyA_ALL_polyA_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
dds_AML_polyA_ALL_polyA_updated
```

    class: DESeqDataSet 
    dim: 25002 40 
    metadata(1): version
    assays(1): counts
    rownames(25002): ENSG00000000003.14 ENSG00000000419.12 ...
      ENSG00000282772.1 ENSG00000282798.1
    rowData names(0):
    colnames(40): THR24_1667_S01 THR24_2131_S01 ... TCGA-AB-3000-03
      TCGA-AB-2921-03
    colData names(1): disease_and_prep

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
AML_polyA_ALL_polyA_updated_genes_to_keep <- rowSums(counts(dds_AML_polyA_ALL_polyA_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_AML_pA_ALL_pA_updated <- dds_AML_polyA_ALL_polyA_updated[which(AML_polyA_ALL_polyA_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_AML_pA_ALL_pA_updated_object <- DESeq(dds_AML_pA_ALL_pA_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 1120 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
result_dds_AML_pA_ALL_pA_updated <- results(dds_AML_pA_ALL_pA_updated_object)
# summarize results table
summary(result_dds_AML_pA_ALL_pA_updated, alpha = 0.05)
```


    out of 24897 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 8783, 35%
    LFC < 0 (down)     : 7786, 31%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.004%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# convert results to dataframe
result_dds_AML_pA_ALL_pA_updated_df <- result_dds_AML_pA_ALL_pA_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
AML_polyA_ALL_polyA_LFC1great_padj <- result_dds_AML_pA_ALL_pA_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
AML_polyA_ALL_polyA_LFC1great_padj
```

    [1] 5586

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
AML_polyA_ALL_polyA_LFC1less_padj <- result_dds_AML_pA_ALL_pA_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
AML_polyA_ALL_polyA_LFC1less_padj
```

    [1] 3878

``` r
hugo_results_AML_polyA_ALL_polyA_updated <- result_dds_AML_pA_ALL_pA_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

# To check whether we lost any DEGs during the HugoID conversion,
# I repeated the above chunk to count the DEGs after converting to HugoIDs
# the numbers were the same

write_tsv(hugo_results_AML_polyA_ALL_polyA_updated, path = "../../output_data/ALL_AML/hugo_results_AML_polyA_ALL_polyA_updated.tsv.gz")
```

    Warning: The `path` argument of `write_tsv()` is deprecated as of readr 1.4.0.
    ℹ Please use the `file` argument instead.

### **AML riboD vs ALL riboD**

``` r
# metadata
AML_riboD_ALL_riboD_meta <- filtered_ALL_AML_list %>%
  filter(disease_and_prep == "AML_riboD" | disease_and_prep == "ALL_riboD")

# expression
AML_riboD_ALL_riboD_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", AML_riboD_ALL_riboD_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_AML_riboD_ALL_riboD_updated_expr <- AML_riboD_ALL_riboD_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_AML_riboD_ALL_riboD_meta <- AML_riboD_ALL_riboD_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_AML_riboD_ALL_riboD_meta) %in% colnames(fix_AML_riboD_ALL_riboD_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_AML_riboD_ALL_riboD_meta) == colnames(fix_AML_riboD_ALL_riboD_updated_expr))
```

    [1] TRUE

``` r
# format data
round_AML_riboD_ALL_riboD_updated_expr = round(fix_AML_riboD_ALL_riboD_updated_expr)

AML_riboD_ALL_riboD_updated_expr_matrix <- as.matrix(round_AML_riboD_ALL_riboD_updated_expr)

fix_AML_riboD_ALL_riboD_meta$disease_and_prep <- as.factor(fix_AML_riboD_ALL_riboD_meta$disease_and_prep)
```

``` r
# set reference levels
# we want to make sure the references are consistent across the comparisons
levels(fix_AML_riboD_ALL_riboD_meta$disease_and_prep)
```

    [1] "ALL_riboD" "AML_riboD"

``` r
# set up DESeq function
dds_AML_riboD_ALL_riboD_updated <- DESeqDataSetFromMatrix(countData = AML_riboD_ALL_riboD_updated_expr_matrix,
                              colData = fix_AML_riboD_ALL_riboD_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
AML_riboD_ALL_riboD_updated_genes_to_keep <- rowSums(counts(dds_AML_riboD_ALL_riboD_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_AML_riboD_ALL_riboD_updated <- dds_AML_riboD_ALL_riboD_updated[which(AML_riboD_ALL_riboD_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_AML_riboD_ALL_riboD_updated_object <- DESeq(dds_AML_riboD_ALL_riboD_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 832 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
results_AML_riboD_ALL_riboD_updated <- results(dds_AML_riboD_ALL_riboD_updated_object)

# summarize results table
summary(results_AML_riboD_ALL_riboD_updated, alpha = 0.05)
```


    out of 24990 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 6650, 27%
    LFC < 0 (down)     : 6547, 26%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.004%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# first we convert to df
results_AML_riboD_ALL_riboD_updated_df <- results_AML_riboD_ALL_riboD_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
AML_riboD_ALL_riboD_LFC1great_padj <- results_AML_riboD_ALL_riboD_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
AML_riboD_ALL_riboD_LFC1great_padj
```

    [1] 3599

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
AML_riboD_ALL_riboD_LFC1less_padj <- results_AML_riboD_ALL_riboD_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
AML_riboD_ALL_riboD_LFC1less_padj
```

    [1] 3116

``` r
hugo_results_AML_riboD_ALL_riboD_updated <- results_AML_riboD_ALL_riboD_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_AML_riboD_ALL_riboD_updated, path = "../../output_data/ALL_AML/hugo_results_AML_riboD_ALL_riboD_updated.tsv.gz")
```

### **AML polyA vs ALL riboD**

``` r
# metadata
AML_polyA_ALL_riboD_meta <- filtered_ALL_AML_list %>%
  filter(disease_and_prep == "AML_polyA" | disease_and_prep == "ALL_riboD")
# the pipe | means we're taking all samples called ALL_riboD OR ALL_polyA in column disease_and_prep

# expression
AML_polyA_ALL_riboD_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", AML_polyA_ALL_riboD_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_AML_polyA_ALL_riboD_updated_expr <- AML_polyA_ALL_riboD_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_AML_polyA_ALL_riboD_meta <- AML_polyA_ALL_riboD_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_AML_polyA_ALL_riboD_meta) %in% colnames(fix_AML_polyA_ALL_riboD_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_AML_polyA_ALL_riboD_meta) == colnames(fix_AML_polyA_ALL_riboD_updated_expr))
```

    [1] TRUE

``` r
# format data
round_AML_polyA_ALL_riboD_updated_expr = round(fix_AML_polyA_ALL_riboD_updated_expr)

AML_polyA_ALL_riboD_updated_expr_matrix <- as.matrix(round_AML_polyA_ALL_riboD_updated_expr)

fix_AML_polyA_ALL_riboD_meta$disease_and_prep <- as.factor(fix_AML_polyA_ALL_riboD_meta$disease_and_prep)
```

``` r
# check reference levels
levels(fix_AML_polyA_ALL_riboD_meta$disease_and_prep)
```

    [1] "ALL_riboD" "AML_polyA"

``` r
# set up DESeq function
dds_AML_polyA_ALL_riboD_updated <- DESeqDataSetFromMatrix(countData = AML_polyA_ALL_riboD_updated_expr_matrix,
                              colData = fix_AML_polyA_ALL_riboD_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
dds_AML_polyA_ALL_riboD_updated
```

    class: DESeqDataSet 
    dim: 25002 40 
    metadata(1): version
    assays(1): counts
    rownames(25002): ENSG00000000003.14 ENSG00000000419.12 ...
      ENSG00000282772.1 ENSG00000282798.1
    rowData names(0):
    colnames(40): THR24_4203_S01 THR24_3471_S01 ... TCGA-AB-3000-03
      TCGA-AB-2921-03
    colData names(1): disease_and_prep

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
AML_polyA_ALL_riboD_updated_genes_to_keep <- rowSums(counts(dds_AML_polyA_ALL_riboD_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_AML_polyA_ALL_riboD_updated <- dds_AML_polyA_ALL_riboD_updated[which(AML_polyA_ALL_riboD_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_AML_polyA_ALL_riboD_updated_object <- DESeq(dds_AML_polyA_ALL_riboD_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 702 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
results_dds_AML_polyA_ALL_riboD_updated <- results(dds_AML_polyA_ALL_riboD_updated_object)
# summarize results table
summary(results_dds_AML_polyA_ALL_riboD_updated, alpha = 0.05)
```


    out of 25000 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 8596, 34%
    LFC < 0 (down)     : 9952, 40%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.004%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# convert results to dataframe
results_dds_AML_polyA_ALL_riboD_updated_df <- results_dds_AML_polyA_ALL_riboD_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
AML_polyA_ALL_riboD_LFC1great_padj <- results_dds_AML_polyA_ALL_riboD_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
AML_polyA_ALL_riboD_LFC1great_padj
```

    [1] 5140

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
AML_polyA_ALL_riboD_LFC1less_padj <- results_dds_AML_polyA_ALL_riboD_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
AML_polyA_ALL_riboD_LFC1less_padj
```

    [1] 6346

``` r
hugo_results_AML_polyA_ALL_riboD_updated <- results_dds_AML_polyA_ALL_riboD_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_AML_polyA_ALL_riboD_updated, path = "../../output_data/ALL_AML/hugo_results_AML_polyA_ALL_riboD_updated.tsv.gz")
```

### **AML riboD vs ALL polyA**

``` r
# metadata
AML_riboD_ALL_polyA_meta <- filtered_ALL_AML_list %>%
  filter(disease_and_prep == "AML_riboD" | disease_and_prep == "ALL_polyA")

# expression
AML_riboD_ALL_polyA_updated_expr <- updated_expression %>%
  select(all_of(c("Gene", AML_riboD_ALL_polyA_meta$term)))
```

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_AML_riboD_ALL_polyA_updated_expr <- AML_riboD_ALL_polyA_updated_expr %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names for clin_merged_STS
fix_AML_riboD_ALL_polyA_meta <- AML_riboD_ALL_polyA_meta %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if row names are same as column names
all(rownames(fix_AML_riboD_ALL_polyA_meta) %in% colnames(fix_AML_riboD_ALL_polyA_updated_expr))
```

    [1] TRUE

``` r
# check in order is the same
all(rownames(fix_AML_riboD_ALL_polyA_meta) == colnames(fix_AML_riboD_ALL_polyA_updated_expr))
```

    [1] TRUE

``` r
# format data
round_AML_riboD_ALL_polyA_updated_expr = round(fix_AML_riboD_ALL_polyA_updated_expr)

AML_riboD_ALL_polyA_updated_expr_matrix <- as.matrix(round_AML_riboD_ALL_polyA_updated_expr)

fix_AML_riboD_ALL_polyA_meta$disease_and_prep <- as.factor(fix_AML_riboD_ALL_polyA_meta$disease_and_prep)
```

``` r
# check reference levels
levels(fix_AML_riboD_ALL_polyA_meta$disease_and_prep)
```

    [1] "ALL_polyA" "AML_riboD"

``` r
# set up DESeq function
dds_AML_riboD_ALL_polyA_updated <- DESeqDataSetFromMatrix(countData = AML_riboD_ALL_polyA_updated_expr_matrix,
                              colData = fix_AML_riboD_ALL_polyA_meta,
                              design = ~disease_and_prep)
```

    converting counts to integer mode

``` r
# create a vector of TRUE and FALSE values where
# TRUE corresponds to genes with counts of at least 10 
AML_riboD_ALL_polyA_updated_genes_to_keep <- rowSums(counts(dds_AML_riboD_ALL_polyA_updated)) >= 10

# use which() to prevent any NAs sneaking through
dds_AML_riboD_ALL_polyA_updated <- dds_AML_riboD_ALL_polyA_updated[which(AML_riboD_ALL_polyA_updated_genes_to_keep), ]
```

``` r
# Run DESeq2
dds_AML_riboD_ALL_polyA_updated_object <- DESeq(dds_AML_riboD_ALL_polyA_updated)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 1169 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
# create results table
results_AML_riboD_ALL_polyA_updated <- results(dds_AML_riboD_ALL_polyA_updated_object)

# summarize results table
summary(results_AML_riboD_ALL_polyA_updated, alpha = 0.05)
```


    out of 25000 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 10213, 41%
    LFC < 0 (down)     : 7472, 30%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

``` r
# first we convert to df
results_AML_riboD_ALL_polyA_updated_df <- results_AML_riboD_ALL_polyA_updated %>%
  #convert to a data frame
  as.data.frame() %>%
  #the gene ids were stored as row names -- let's them a column
  tibble::rownames_to_column(var = "Gene")
```

``` r
# how many genes have a Log2 fold change greater or equal to 1 and have an adjusted p-value less than 0.05?
AML_riboD_ALL_polyA_LFC1great_padj <- results_AML_riboD_ALL_polyA_updated_df %>%
  filter(log2FoldChange >= 1 & padj < 0.05) %>%
  nrow()
AML_riboD_ALL_polyA_LFC1great_padj
```

    [1] 7252

``` r
# how many genes have a Log2 fold change less than or equal to -1 and have an adjusted p-value less than 0.05?
AML_riboD_ALL_polyA_LFC1less_padj <- results_AML_riboD_ALL_polyA_updated_df %>%
  filter(log2FoldChange <= -1 & padj < 0.05) %>%
  nrow()
AML_riboD_ALL_polyA_LFC1less_padj
```

    [1] 4115

``` r
hugo_results_AML_riboD_ALL_polyA_updated <- results_AML_riboD_ALL_polyA_updated_df %>%
  left_join(gene_names, by=c("Gene"="EnsGeneID")) %>%
  relocate(HugoID) # moving HugoID column to the front

write_tsv(hugo_results_AML_riboD_ALL_polyA_updated, path = "../../output_data/ALL_AML/hugo_results_AML_riboD_ALL_polyA_updated.tsv.gz")
```

#### Session Info

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
     date     2026-07-02
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

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/Fig_2/ALL_AML/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Users/maryke/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
