# Fig_1C


## Fig 1C

Stacked bar chart showing the proportion of unique genes expressed at
different log2(TPM+1) levels.

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ✔ purrr     1.2.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# expression files
SS_aRMS_log2tpm1 <- read_tsv("../input_data/log2tpm_ensembl_SS_aRMS.tsv.gz")
```

    Rows: 60498 Columns: 61
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (60): THR51_4556_S01, THR51_4558_S01, THR51_4554_S01, THR24_3992_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_NB_log2tpm1 <- read_tsv("../input_data/log2tpm_ensembl_WT_NB.tsv.gz")
```

    Rows: 60498 Columns: 61
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (60): THR24_3218_S01, THR24_4284_S01, THR24_4194_S01, THR24_4369_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_AML_log2tpm1 <- read_tsv("../input_data/log2tpm_ensembl_ALL_AML.tsv.gz")
```

    Rows: 60498 Columns: 81
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (80): THR24_4203_S01, THR24_3471_S01, THR24_4235_S01, THR24_3688_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# sample ID files, from correlation_analysis
filtered_SS_aRMS_list <- read_tsv("../input_data/filtered_SS_aRMS_list.tsv")
```

    Rows: 60 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
filtered_WT_NB_list <- read_tsv("../input_data/filtered_WT_NB_list.tsv")
```

    Rows: 60 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
filtered_ALL_AML_list <- read_tsv("../input_data/filtered_ALL_AML_list.tsv")
```

    Rows: 80 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Since "Gene" is a column name, it won't match the row names
# use values from the first column ("Gene") as row names
fix_SS_aRMS_log2tpm1 <- SS_aRMS_log2tpm1 %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

fix_WT_NB_log2tpm1 <- WT_NB_log2tpm1 %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

fix_ALL_AML_log2tpm1 <- ALL_AML_log2tpm1 %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")

# use values from first column ("th_dataset_id") as row names
fix_filtered_SS_aRMS_list <- filtered_SS_aRMS_list %>%
  remove_rownames %>%
  column_to_rownames(var = "term")

fix_filtered_WT_NB_list <- filtered_WT_NB_list %>%
  remove_rownames %>%
  column_to_rownames(var = "term")

fix_filtered_ALL_AML_list <- filtered_ALL_AML_list %>%
  remove_rownames %>%
  column_to_rownames(var = "term")
```

``` r
# check if samples are the same in expression and ID files
all(rownames(fix_filtered_SS_aRMS_list) %in% colnames(fix_SS_aRMS_log2tpm1))
```

    [1] TRUE

``` r
all(rownames(fix_filtered_WT_NB_list) %in% colnames(fix_WT_NB_log2tpm1))
```

    [1] FALSE

``` r
all(rownames(fix_filtered_ALL_AML_list) %in% colnames(fix_ALL_AML_log2tpm1))
```

    [1] FALSE

``` r
# Returns row names not found in column names
setdiff(rownames(fix_filtered_WT_NB_list), colnames(fix_WT_NB_log2tpm1))
```

    [1] "TARGET-30-PASPBZ-01" "TARGET-30-PATYIL-02" "TARGET-30-PARHAM-01"

``` r
# Returns column names not found in row names
setdiff(colnames(fix_WT_NB_log2tpm1), rownames(fix_filtered_WT_NB_list))
```

    [1] "TARGET-30-PATGLU-01" "TARGET-30-PARACS-01" "TARGET-30-PAPLSD-01"

``` r
# Returns row names not found in column names
setdiff(rownames(fix_filtered_ALL_AML_list), colnames(fix_ALL_AML_log2tpm1))
```

    [1] "THR24_2009_S01"  "THR24_1704_S01"  "THR24_1670_S01"  "THR24_1887_S01" 
    [5] "THR24_4159_S01"  "THR24_3433_S01"  "TCGA-AB-2921-03" "THR24_3307_S01" 

``` r
# Returns column names not found in row names
setdiff(colnames(fix_ALL_AML_log2tpm1), rownames(fix_filtered_ALL_AML_list))
```

    [1] "THR24_3668_S01"  "THR24_3660_S01"  "THR24_2112_S01"  "THR24_1581_S01" 
    [5] "THR24_2132_S01"  "THR24_1570_S01"  "TCGA-AB-2930-03" "THR24_4297_S01" 
