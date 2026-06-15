# Fig1G


## Fig 1G

Scatterplot showing median expression per gene.

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
SS_aRMS_log2tpm1 <- read_tsv("../../input_data/SS_aRMS_log2TPM1_ensembl_TEBA.tsv.gz")
```

    Rows: 60498 Columns: 61
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (60): THR51_4556_S01, THR51_4558_S01, THR51_4554_S01, THR24_3992_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_NB_log2tpm1 <- read_tsv("../../input_data/WT_NB_log2TPM1_ensembl_TEBA.tsv.gz")
```

    Rows: 60498 Columns: 61
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (60): THR24_3218_S01, THR24_4194_S01, THR24_4284_S01, THR24_4369_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_AML_log2tpm1 <- read_tsv("../../input_data/ALL_AML_log2TPM1_ensembl_TEBA.tsv.gz")
```

    Rows: 60498 Columns: 81
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (80): THR24_1667_S01, THR24_2131_S01, THR24_2119_S01, THR24_1921_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# sample ID files, from correlation_analysis
filtered_SS_aRMS_list <- read_tsv("../../input_data/filtered_SS_aRMS_list.tsv")
```

    Rows: 60 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
filtered_WT_NB_list <- read_tsv("../../input_data/filtered_WT_NB_list.tsv")
```

    Rows: 60 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
filtered_ALL_AML_list <- read_tsv("../../input_data/filtered_ALL_AML_list.tsv")
```

    Rows: 80 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# to convert EnsemblIDs to HugoIDs
gene_names <- read.table("../../input_data/EnsGeneID_Hugo_Observed_Conversions.txt",
header = TRUE, sep = "\t", stringsAsFactors = FALSE )
```

``` r
# SS polyA
SS_polyA_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("SS_polyA")) %>% 
  select(term, disease_and_prep)

SS_polyA_expr <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(SS_polyA_list$term))

# SS riboD
SS_riboD_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("SS_riboD")) %>% 
  select(term, disease_and_prep)

SS_riboD_expr <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(SS_riboD_list$term))

# aRMS polyA
aRMS_polyA_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("aRMS_polyA")) %>% 
  select(term, disease_and_prep)

aRMS_polyA_expr <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(aRMS_polyA_list$term))

# aRMS riboD
aRMS_riboD_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("aRMS_riboD")) %>% 
  select(term, disease_and_prep)

aRMS_riboD_expr <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(aRMS_riboD_list$term))

# WT polyA
WT_polyA_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("WT_polyA")) %>% 
  select(term, disease_and_prep)

WT_polyA_expr <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(WT_polyA_list$term))

# WT riboD 
WT_riboD_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("WT_riboD")) %>% 
  select(term, disease_and_prep)

WT_riboD_expr <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(WT_riboD_list$term))

# NB polyA
NB_polyA_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("NB_polyA")) %>% 
  select(term, disease_and_prep)

NB_polyA_expr <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(NB_polyA_list$term))

# NB riboD
NB_riboD_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("NB_riboD")) %>% 
  select(term, disease_and_prep)

NB_riboD_expr <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(NB_riboD_list$term))

# ALL polyA
ALL_polyA_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("ALL_polyA")) %>% 
  select(term, disease_and_prep)

ALL_polyA_expr <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(ALL_polyA_list$term))

# ALL riboD
ALL_riboD_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("ALL_riboD")) %>% 
  select(term, disease_and_prep)

ALL_riboD_expr <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(ALL_riboD_list$term))


# AML polyA
AML_polyA_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("AML_polyA")) %>% 
  select(term, disease_and_prep)

AML_polyA_expr <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(AML_polyA_list$term))

# AML riboD
AML_riboD_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("AML_riboD")) %>% 
  select(term, disease_and_prep)

AML_riboD_expr <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(AML_riboD_list$term))
```
