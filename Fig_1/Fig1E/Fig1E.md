# Fig1E


## Fig 1E

Boxplot showing expression level of HIST1H1B, a nonpolyadenylated
histone subunit gene.

Boxplot showing expression of average median ratio gene, calculated as
the median expression level in PolyA samples divided by the median
expression level in RiboD samples

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

Convert EnsemblIDs to HugoIDs

``` r
convert_to_hugo <- function(df) {
  df %>%
    left_join(gene_names, by = c("Gene" = "EnsGeneID")) %>%
    mutate(Gene = ifelse(!is.na(HugoID), HugoID, Gene)) %>%
    select(-HugoID)
}
```

Finding median expression genes

``` r
# Find gene medians in polyA SS
SS_polyA_medians <- SS_polyA_expr %>%
  rowwise() %>%
  mutate(SS_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, SS_polyA_median)

# Find gene medians in riboD SS
SS_riboD_medians <- SS_riboD_expr %>%
  rowwise() %>%
  mutate(SS_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, SS_riboD_median)

# Join and compute ratio (polyA median/riboD median)
SS_gene_medians <- SS_polyA_medians %>%
  inner_join(SS_riboD_medians, by = "Gene") %>%
  mutate(SS_median_ratio = SS_polyA_median / SS_riboD_median)

# Filter out genes with no expression
SS_gene_medians_filtered <- SS_gene_medians %>%
  filter(SS_polyA_median > 0 & SS_riboD_median > 0)

# Recompute median ratio value
SS_median_ratio_val <- median(SS_gene_medians_filtered$SS_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
SS_most_average_gene <- SS_gene_medians_filtered %>%
  mutate(dist = abs(SS_median_ratio - SS_median_ratio_val)) %>%
  arrange(dist) %>%
  slice(1)  

SS_most_average_gene_hugo <- SS_most_average_gene %>%
  convert_to_hugo()
SS_most_average_gene_hugo
```

| Gene  | SS_polyA_median | SS_riboD_median | SS_median_ratio |  dist |
|:------|----------------:|----------------:|----------------:|------:|
| MTMR2 |        4.797532 |        2.954196 |        1.623972 | 5e-07 |

``` r
# Find gene medians in polyA aRMS
aRMS_polyA_medians <- aRMS_polyA_expr %>%
  rowwise() %>%
  mutate(aRMS_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, aRMS_polyA_median)

# Find gene medians in riboD aRMS
aRMS_riboD_medians <- aRMS_riboD_expr %>%
  rowwise() %>%
  mutate(aRMS_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, aRMS_riboD_median)

# Join and compute ratio (polyA median/riboD median)
aRMS_gene_medians <- aRMS_polyA_medians %>%
  inner_join(aRMS_riboD_medians, by = "Gene") %>%
  mutate(aRMS_median_ratio = aRMS_polyA_median / aRMS_riboD_median)

# Filter out genes with no or near-zero expression
aRMS_gene_medians_filtered <- aRMS_gene_medians %>%
  filter(aRMS_polyA_median > 0 & aRMS_riboD_median > 0)

# Recompute median ratio value
aRMS_median_ratio_val <- median(aRMS_gene_medians_filtered$aRMS_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
aRMS_most_average_gene <- aRMS_gene_medians_filtered %>%
  mutate(dist = abs(aRMS_median_ratio - aRMS_median_ratio_val)) %>%  # computes how far each gene’s ratio is from overall median
  arrange(dist) %>%
  slice(1)   # Isolates gene closest to ave median ration

aRMS_most_average_gene_hugo <- aRMS_most_average_gene %>%
  convert_to_hugo()
aRMS_most_average_gene_hugo
```

| Gene  | aRMS_polyA_median | aRMS_riboD_median | aRMS_median_ratio | dist |
|:------|------------------:|------------------:|------------------:|-----:|
| TM2D1 |           5.01078 |          3.065228 |          1.634717 |    0 |

``` r
# Find gene medians in polyA WT
WT_polyA_medians <- WT_polyA_expr %>%
  rowwise() %>%
  mutate(WT_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, WT_polyA_median)

# Find gene medians in riboD WT
WT_riboD_medians <- WT_riboD_expr %>%
  rowwise() %>%
  mutate(WT_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, WT_riboD_median)

# Join and compute ratio (polyA median/riboD median)
WT_gene_medians <- WT_polyA_medians %>%
  inner_join(WT_riboD_medians, by = "Gene") %>%
  mutate(WT_median_ratio = WT_polyA_median / WT_riboD_median)

# Filter out genes with no or near-zero expression
WT_gene_medians_filtered <- WT_gene_medians %>%
  filter(WT_polyA_median > 0 & WT_riboD_median > 0)

# Recompute median ratio value
WT_median_ratio_val <- median(WT_gene_medians_filtered$WT_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
WT_most_average_gene <- WT_gene_medians_filtered %>%
  mutate(dist = abs(WT_median_ratio - WT_median_ratio_val)) %>%
  arrange(dist) %>%
  slice(1)         # Isolates gene closest to ave median ration

WT_most_average_gene_hugo <- WT_most_average_gene %>%
  convert_to_hugo()
WT_most_average_gene_hugo
```

| Gene  | WT_polyA_median | WT_riboD_median | WT_median_ratio |     dist |
|:------|----------------:|----------------:|----------------:|---------:|
| SYPL1 |        6.028748 |        3.587365 |         1.68055 | 4.04e-05 |

``` r
# Find gene medians in polyA NB
NB_polyA_medians <- NB_polyA_expr %>%
  rowwise() %>%
  mutate(NB_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, NB_polyA_median)

# Find gene medians in riboD NB
NB_riboD_medians <- NB_riboD_expr %>%
  rowwise() %>%
  mutate(NB_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, NB_riboD_median)

# Join and compute ratio (polyA median/riboD median)
NB_gene_medians <- NB_polyA_medians %>%
  inner_join(NB_riboD_medians, by = "Gene") %>%
  mutate(NB_median_ratio = NB_polyA_median / NB_riboD_median)

# Filter out genes with no or near-zero expression
NB_gene_medians_filtered <- NB_gene_medians %>%
  filter(NB_polyA_median > 0 & NB_riboD_median > 0)

# Recompute median ratio value
NB_median_ratio_val <- median(NB_gene_medians_filtered$NB_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
NB_most_average_gene <- NB_gene_medians_filtered %>%
  mutate(dist = abs(NB_median_ratio - NB_median_ratio_val)) %>%  # computes how far each gene’s ratio is from overall median
  arrange(dist) %>%
  slice(1)         # Isolates gene closest to ave median ration

NB_most_average_gene_hugo <- NB_most_average_gene %>%
  convert_to_hugo()
NB_most_average_gene_hugo
```

| Gene   | NB_polyA_median | NB_riboD_median | NB_median_ratio |     dist |
|:-------|----------------:|----------------:|----------------:|---------:|
| CCDC64 |        5.667639 |        3.228049 |        1.755747 | 5.88e-05 |

``` r
# Find gene medians in polyA ALL
ALL_polyA_medians <- ALL_polyA_expr %>%
  rowwise() %>%
  mutate(ALL_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, ALL_polyA_median)

# Find gene medians in riboD ALL
ALL_riboD_medians <- ALL_riboD_expr %>%
  rowwise() %>%
  mutate(ALL_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, ALL_riboD_median)

# Join and compute ratio (polyA median/riboD median)
ALL_gene_medians <- ALL_polyA_medians %>%
  inner_join(ALL_riboD_medians, by = "Gene") %>%
  mutate(ALL_median_ratio = ALL_polyA_median / ALL_riboD_median)

# Filter out genes with no or near-zero expression
ALL_gene_medians_filtered <- ALL_gene_medians %>%
  filter(ALL_polyA_median > 0 & ALL_riboD_median > 0)

# Recompute median ratio value
ALL_median_ratio_val <- median(ALL_gene_medians_filtered$ALL_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
ALL_most_average_gene <- ALL_gene_medians_filtered %>%
  mutate(dist = abs(ALL_median_ratio - ALL_median_ratio_val)) %>%  # computes how far each gene’s ratio is from overall median
  arrange(dist) %>%
  slice(1)         # Isolates gene closest to ave median ration

ALL_most_average_gene_hugo <- ALL_most_average_gene %>%
  convert_to_hugo()
ALL_most_average_gene_hugo
```

| Gene | ALL_polyA_median | ALL_riboD_median | ALL_median_ratio |    dist |
|:-----|-----------------:|-----------------:|-----------------:|--------:|
| XPO5 |          5.28429 |         2.122672 |         2.489452 | 1.3e-05 |

``` r
# Find gene medians in polyA AML
AML_polyA_medians <- AML_polyA_expr %>%
  rowwise() %>%
  mutate(AML_polyA_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, AML_polyA_median)

# Find gene medians in riboD AML
AML_riboD_medians <- AML_riboD_expr %>%
  rowwise() %>%
  mutate(AML_riboD_median = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Gene, AML_riboD_median)

# Join and compute ratio (polyA median/riboD median)
AML_gene_medians <- AML_polyA_medians %>%
  inner_join(AML_riboD_medians, by = "Gene") %>%
  mutate(AML_median_ratio = AML_polyA_median / AML_riboD_median)

# Filter out genes with no or near-zero expression
AML_gene_medians_filtered <- AML_gene_medians %>%
  filter(AML_polyA_median > 1 & AML_riboD_median > 0.9)

# Recompute median ratio value
AML_median_ratio_val <- median(AML_gene_medians_filtered$AML_median_ratio, na.rm = TRUE)

# Find gene closest to that median ratio
AML_most_average_gene <- AML_gene_medians_filtered %>%
  mutate(dist = abs(AML_median_ratio - AML_median_ratio_val)) %>%  # computes how far each gene’s ratio is from overall median
  arrange(dist) %>%
  slice(1)         # Isolates gene closest to ave median ration

AML_most_average_gene_hugo <- AML_most_average_gene %>%
  convert_to_hugo()
AML_most_average_gene_hugo
```

| Gene | AML_polyA_median | AML_riboD_median | AML_median_ratio |     dist |
|:-----|-----------------:|-----------------:|-----------------:|---------:|
| EML3 |         5.895521 |         3.329959 |         1.770449 | 7.04e-05 |

Converting expression to hugoIDs

``` r
# SS
SS_polyA_expression_hugo <- SS_polyA_expr %>%
  convert_to_hugo()

SS_riboD_expression_hugo <- SS_riboD_expr %>%
  convert_to_hugo()

SS_gene_medians_hugo <- SS_gene_medians %>%
  convert_to_hugo()

# aRMS
aRMS_polyA_expression_hugo <- aRMS_polyA_expr %>%
  convert_to_hugo() 

aRMS_riboD_expression_hugo <- aRMS_riboD_expr %>%
  convert_to_hugo()

aRMS_gene_median_hugo <- aRMS_gene_medians %>%
  convert_to_hugo()

# WT
WT_polyA_expression_hugo <- WT_polyA_expr %>%
  convert_to_hugo()

WT_riboD_expression_hugo <- WT_riboD_expr %>%
  convert_to_hugo()

WT_gene_medians_hugo <- WT_gene_medians %>%
  convert_to_hugo()

# NB
NB_polyA_expression_hugo <- NB_polyA_expr %>%
  convert_to_hugo()

NB_riboD_expression_hugo <- NB_riboD_expr %>%
  convert_to_hugo()

NB_gene_medians_hugo <- NB_gene_medians %>%
  convert_to_hugo()

# ALL
ALL_polyA_expression_hugo <- ALL_polyA_expr %>%
  convert_to_hugo()

ALL_riboD_expression_hugo <- ALL_riboD_expr %>%
  convert_to_hugo()

ALL_gene_medians_hugo <- ALL_gene_medians %>%
  convert_to_hugo()

# AML
AML_polyA_expression_hugo <- AML_polyA_expr %>%
  convert_to_hugo()

AML_riboD_expression_hugo <- AML_riboD_expr %>%
  convert_to_hugo()

AML_gene_medians_hugo <- AML_gene_medians %>%
  convert_to_hugo()
```

Filtering for HIST1H1B, a non-polyadenylated gene

``` r
H1B_SS_polyA_expression <- SS_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_SS_riboD_expression <- SS_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_aRMS_polyA_expression <- aRMS_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B")

H1B_aRMS_riboD_expression <- aRMS_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_WT_polyA_expression <- WT_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B")

H1B_WT_riboD_expression <- WT_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_NB_polyA_expression <- NB_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B")

H1B_NB_riboD_expression <- NB_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_ALL_polyA_expression <- ALL_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B")

H1B_ALL_riboD_expression <- ALL_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 

H1B_AML_polyA_expression <- AML_polyA_expression_hugo %>%
   filter(Gene == "HIST1H1B")

H1B_AML_riboD_expression <- AML_riboD_expression_hugo %>%
   filter(Gene == "HIST1H1B") 
```

Filtering for most averagely expressed gene

``` r
#SS Most Average Expressed gene
avg_SS_polyA_expression <- SS_polyA_expression_hugo %>%
   filter(Gene == SS_most_average_gene_hugo$Gene) 

avg_SS_riboD_expression <- SS_riboD_expression_hugo %>%
   filter(Gene == SS_most_average_gene_hugo$Gene)

#aRMS Most Average Expressed gene
avg_aRMS_riboD_expression <- aRMS_riboD_expression_hugo %>%
   filter(Gene == aRMS_most_average_gene_hugo$Gene)

avg_aRMS_polyA_expression <- aRMS_polyA_expression_hugo %>%
   filter(Gene == aRMS_most_average_gene_hugo$Gene)

#WT Most Average Expressed gene 
avg_WT_riboD_expression <- WT_riboD_expression_hugo %>%
   filter(Gene == WT_most_average_gene_hugo$Gene)

avg_WT_polyA_expression <- WT_polyA_expression_hugo %>%
   filter(Gene == WT_most_average_gene_hugo$Gene)

#NB Most Average Expressed gene
avg_NB_riboD_expression <- NB_riboD_expression_hugo %>%
   filter(Gene == NB_most_average_gene_hugo$Gene)

avg_NB_polyA_expression <- NB_polyA_expression_hugo %>%
   filter(Gene == NB_most_average_gene_hugo$Gene)

#ALL Most Average Expressed gene
avg_ALL_riboD_expression <- ALL_riboD_expression_hugo %>%
   filter(Gene == ALL_most_average_gene_hugo$Gene)

avg_ALL_polyA_expression <- ALL_polyA_expression_hugo %>%
   filter(Gene == ALL_most_average_gene_hugo$Gene)

#AML Most Average Expressed gene
avg_AML_riboD_expression <- AML_riboD_expression_hugo %>%
   filter(Gene == AML_most_average_gene_hugo$Gene)

avg_AML_polyA_expression <- AML_polyA_expression_hugo %>%
   filter(Gene == AML_most_average_gene_hugo$Gene)
```

Reformating HIST1H1B data

``` r
H1B_SS_polyA_long <- H1B_SS_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "SS", Compendia = "polyA", Gene = "HIST1H1B")

H1B_SS_riboD_long <- H1B_SS_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "SS", Compendia = "riboD", Gene = "HIST1H1B")

H1B_aRMS_polyA_long <- H1B_aRMS_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "aRMS", Compendia = "polyA", Gene = "HIST1H1B")

H1B_aRMS_riboD_long <- H1B_aRMS_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "aRMS", Compendia = "riboD", Gene = "HIST1H1B")

H1B_WT_polyA_long <- H1B_WT_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "WT", Compendia = "polyA", Gene = "HIST1H1B")

H1B_WT_riboD_long <- H1B_WT_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "WT", Compendia = "riboD", Gene = "HIST1H1B")

H1B_NB_polyA_long <- H1B_NB_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "NB", Compendia = "polyA", Gene = "HIST1H1B")

H1B_NB_riboD_long <- H1B_NB_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "NB", Compendia = "riboD", Gene = "HIST1H1B")

H1B_ALL_polyA_long <- H1B_ALL_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "ALL", Compendia = "polyA", Gene = "HIST1H1B")

H1B_ALL_riboD_long <- H1B_ALL_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "ALL", Compendia = "riboD", Gene = "HIST1H1B")

H1B_AML_polyA_long <- H1B_AML_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "AML", Compendia = "polyA", Gene = "HIST1H1B")

H1B_AML_riboD_long <- H1B_AML_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "AML", Compendia = "riboD", Gene = "HIST1H1B")
```

Reformatting avg gene

``` r
# SS polyA 
avg_SS_polyA_long <- avg_SS_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "SS", Compendia = "polyA", Gene = Gene)

# SS riboD
avg_SS_riboD_long <- avg_SS_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "SS", Compendia = "riboD", Gene = Gene)

# aRMS polyA
avg_aRMS_polyA_long <- avg_aRMS_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "aRMS", Compendia = "polyA", Gene = Gene)

# aRMS riboD
avg_aRMS_riboD_long <- avg_aRMS_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "aRMS", Compendia = "riboD", Gene = Gene)

#NB polyA
avg_NB_polyA_long <- avg_NB_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "NB", Compendia = "polyA", Gene = Gene)

#NB riboD
avg_NB_riboD_long <- avg_NB_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "NB", Compendia = "riboD", Gene = Gene)

#WT polyA
avg_WT_polyA_long <- avg_WT_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "WT", Compendia = "polyA", Gene = Gene)

#WT riboD
avg_WT_riboD_long <- avg_WT_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "WT", Compendia = "riboD", Gene = Gene)

#AML polyA
avg_AML_polyA_long <- avg_AML_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "AML", Compendia = "polyA", Gene = Gene)

#AML riboD
avg_AML_riboD_long <- avg_AML_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "AML", Compendia = "riboD", Gene = Gene)

#ALL polyA
avg_ALL_polyA_long <- avg_ALL_polyA_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "ALL", Compendia = "polyA", Gene = Gene)

#ALL riboD
avg_ALL_riboD_long <- avg_ALL_riboD_expression %>%
  pivot_longer(-Gene, names_to = "TH_ID", values_to = "Expression") %>%
  mutate(Disease = "ALL", Compendia = "riboD", Gene = Gene)
```

Combining dataframes

``` r
combined_long <- bind_rows(

  # -------------------------
  # HIST1H1B — SS, aRMS, WT, NB, ALL, AML
  # -------------------------
  H1B_SS_polyA_long %>%
    mutate(Disease = "SS"),
  H1B_SS_riboD_long %>%
    mutate(Disease = "SS"),
  
  H1B_aRMS_polyA_long %>%
    mutate(Disease = "aRMS"),
  H1B_aRMS_riboD_long %>%
    mutate(Disease = "aRMS"),
  
  H1B_WT_polyA_long %>%
    mutate(Disease = "WT"),
  H1B_WT_riboD_long %>%
    mutate(Disease = "WT"),
  
  H1B_NB_polyA_long %>%
    mutate(Disease = "NB"),
  H1B_NB_riboD_long %>%
    mutate(Disease = "NB"),
  
  H1B_ALL_polyA_long %>%
    mutate(Disease = "ALL"),
  H1B_ALL_riboD_long %>%
    mutate(Disease = "ALL"),
  
  H1B_AML_polyA_long %>%
    mutate(Disease = "AML"),
  H1B_AML_riboD_long %>%
    mutate(Disease = "AML"),
   # -------------------------
  # Median-ratio genes (SS, aRMS, WT, NB, ALL, AML)
  # -------------------------
  avg_SS_polyA_long %>%
    mutate(Disease = "SS"),
  avg_SS_riboD_long %>%
    mutate(Disease = "SS"),

  avg_aRMS_polyA_long %>%
    mutate(Disease = "aRMS"),
  avg_aRMS_riboD_long %>%
    mutate(Disease = "aRMS"),

  avg_WT_polyA_long %>%
    mutate(Disease = "WT"),
  avg_WT_riboD_long %>%
    mutate(Disease = "WT"),

  avg_NB_polyA_long %>%
    mutate(Disease = "NB"),
  avg_NB_riboD_long %>%
    mutate(Disease = "NB"),
  
  avg_ALL_polyA_long %>%
    mutate(Disease = "ALL"),
  avg_ALL_riboD_long %>%
    mutate(Disease = "ALL"),
  
  avg_AML_polyA_long %>%
    mutate(Disease = "AML"),
  avg_AML_riboD_long %>%
    mutate(Disease = "AML"),
)
```

Custom color theme

``` r
color_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      plot.title = element_text(face = "bold", size = base_size * 1.1, hjust = 0.5),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 16), 
      axis.title.x = element_text(angle = 0, hjust = 0.5, size = 17), 
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.margin = margin(10, 10, 10, 10)
    )
}

scale_fill_compendia <- function() {
  scale_fill_manual(values = c(
    "riboD" = "#E69F00",  # yellow
    "polyA" = "#0072B2"   # blue
  ))
}
scale_color_compendia <- function() {
  scale_color_manual(values = c( 
    "riboD" = "#E69F00", 
    "polyA" = "#0072B2" ))
  }
```

Merged plot

``` r
# Map of median-ratio genes
median_gene_map <- list(
  SS      = SS_most_average_gene_hugo$Gene,
  aRMS    = aRMS_most_average_gene_hugo$Gene,
  WT      = WT_most_average_gene_hugo$Gene,
  NB      = NB_most_average_gene_hugo$Gene,
  ALL     = ALL_most_average_gene_hugo$Gene,
  AML     = AML_most_average_gene_hugo$Gene
)

diseases <- unique(combined_long$Disease)
```

Generate merged plot for each disease

``` r
merged_plots <- map(diseases, function(d) {

  median_gene <- median_gene_map[[d]]

  # HIST1H1B expression
  df_h1b <- combined_long %>%
    filter(Disease == d, Gene == "HIST1H1B") %>%
    mutate(PlotType = "HIST1H1B (PolyA-)")

  # Median-ratio gene expression
  df_median <- combined_long %>%
    filter(Disease == d, Gene == median_gene) %>%
    mutate(PlotType = "Median Ratio Gene (PolyA+)")

  # Combine
  df_combined <- bind_rows(df_h1b, df_median)

  ggplot(df_combined, aes(x = Gene, y = Expression, fill = Compendia)) +
    geom_boxplot(
      alpha = 0.7,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      linewidth = 0.4
    ) +
    geom_jitter(
      aes(color = Compendia),
      position = position_jitterdodge(jitter.width = 0.2,
                                      dodge.width = 0.8),
      size = 1.2,
      alpha = 0.6
    ) +
    scale_fill_compendia() +
    scale_color_compendia() +
    coord_cartesian(ylim = c(0, 17)) +
    facet_wrap(~ PlotType, ncol = 2, scales = "free_x") +
    labs(
      title = paste("Expression of Representative PolyA+ and PolyA- Gene in", d),
      x = "Gene",
      y = "Expression log2(TPM+1)"
    ) +
    color_theme()
})

names(merged_plots) <- diseases
merged_plots
```

    $SS

![](Fig1E_files/figure-commonmark/Fig1E_nostars-1.png)


    $aRMS

![](Fig1E_files/figure-commonmark/Fig1E_nostars-2.png)


    $WT

![](Fig1E_files/figure-commonmark/Fig1E_nostars-3.png)


    $NB

![](Fig1E_files/figure-commonmark/Fig1E_nostars-4.png)


    $ALL

![](Fig1E_files/figure-commonmark/Fig1E_nostars-5.png)


    $AML

![](Fig1E_files/figure-commonmark/Fig1E_nostars-6.png)

Statistical test for significance

``` r
choose_test <- function(poly, ribo) {
  
  # Shapiro tests
  p_norm_poly  <- shapiro.test(poly)$p.value
  p_norm_ribo  <- shapiro.test(ribo)$p.value
  
  poly_normal <- p_norm_poly  > 0.05
  ribo_normal <- p_norm_ribo  > 0.05
  
  # Both normal → variance test
  if (poly_normal & ribo_normal) {
    
    p_var <- var.test(poly, ribo)$p.value
    equal_var <- p_var > 0.05
    
    # Normal + equal variance → Student t-test
    if (equal_var) {
      return(list(
        p_value = t.test(poly, ribo, var.equal = TRUE)$p.value,
        test_used = "student_t"
      ))
    }
    
    # Normal + unequal variance → Welch t-test
    else {
      return(list(
        p_value = t.test(poly, ribo, var.equal = FALSE)$p.value,
        test_used = "welch_t"
      ))
    }
  }
  
  # Non-normal → Mann–Whitney U test
  return(list(
    p_value = wilcox.test(poly, ribo, exact = FALSE)$p.value,
    test_used = "mann_whitney"
  ))
}
```

Assigning P-value to significance star

``` r
p_to_stars <- function(p) 
{ vapply(p, function(x) 
{ if (is.na(x)) return("NA") 
  if (x < 0.0001) return("****") 
  if (x < 0.001) return("***") 
  if (x < 0.01) return("**") 
  if (x < 0.05) return("*") 
  "ns" 
  }, FUN.VALUE = character(1)) }
```

Applying to data

``` r
compute_significance_table <- function(df) {
  
  # 1. Clean and deduplicate input
  df_clean <- df %>%
    select(Disease, Gene, Compendia, Expression) %>%
    distinct()
  
  # 2. Collapse to one row per Disease × Gene × Compendia
  df_collapsed <- df_clean %>%
    group_by(Disease, Gene, Compendia) %>%
    summarise(Expression = list(Expression), .groups = "drop")
  
  # 3. Pivot to wide: polyA and riboD in separate columns
  df_wide <- df_collapsed %>%
    tidyr::pivot_wider(
      names_from = Compendia,
      values_from = Expression
    )
  
  # 4. Apply statistical test to each row
  df_results <- df_wide %>%
    rowwise() %>%
    mutate(
      # wrap the returned list so mutate() treats it as a single object
      test = list(choose_test(unlist(polyA), unlist(riboD))),
      p_value = test$p_value,
      test_used = test$test_used,
      stars = p_to_stars(p_value)
    ) %>%
    ungroup() %>%
    select(-test)
  
    # Benjamini–Hochberg FDR correction
  df_results2 <- df_results %>%
    group_by(Disease) %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(stars_adj = p_to_stars(padj))
  
  df_results2
}


sig_results <- compute_significance_table(combined_long)

sig_results
```

| Disease | Gene | polyA | riboD | p_value | test_used | stars | padj | stars_adj |
|:---|:---|:---|:---|---:|:---|:---|---:|:---|
| ALL | HIST1H1B | 0.4854268, 0.1763228, 0.0000000, 0.7136958, 0.9030383, 0.5753123, 0.3219281, 3.5668152, 3.1953476, 0.6415460, 3.0959244, 1.3840498, 0.4436067, 0.6690268, 0.7824086, 3.2570106, 1.9068906, 1.0000000 | 7.499368, 7.899659, 7.974415, 6.757157, 6.547974, 7.524189, 6.972922, 7.657497, 6.785681, 8.443441, 7.566130, 8.312611, 7.493855, 7.279750, 7.029232, 7.887099, 7.563463, 7.077884, 6.657640, 6.593802 | 2.00e-07 | mann_whitney | \*\*\*\* | 2.00e-07 | \*\*\*\* |
| ALL | XPO5 | 5.447579, 4.785027, 4.890933, 5.489928, 5.286512, 5.469235, 5.282069, 4.614710, 5.418190, 4.861955, 4.694323, 5.549669, 5.392317, 5.855242, 4.859473, 5.761019, 5.137504, 5.384395, 5.056584, 5.187451 | 1.678072, 2.257011, 2.400538, 1.752749, 2.358959, 2.336283, 2.035624, 2.389567, 2.032101, 2.189034, 2.121015, 2.622930, 1.970854, 2.182692, 2.025029, 1.769772, 2.244887, 2.124328, 1.782409 | 0.00e+00 | student_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| AML | EML3 | 5.421920, 5.820926, 5.648427, 5.876250, 5.752488, 5.950892, 5.330182, 5.914792, 6.077208, 5.600517, 5.498247, 6.350671, 6.094450, 6.376948, 5.727661, 6.108148, 5.389210, 5.928365, 5.952564, 6.134465 | 3.026800, 4.209453, 3.144046, 3.001802, 3.246408, 3.823749, 3.872829, 3.285402, 3.582556, 4.137504, 2.865919, 4.156235, 3.371559, 3.221877, 2.899176, 3.928844, 3.288359, 3.913608, 3.500802, 2.887525 | 0.00e+00 | student_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| AML | HIST1H1B | 0.1244596, 0.4437205, 0.6136316, 0.8156711, 0.0000000, 0.1376322, 0.2631601, 0.4006362, 1.2869577, 0.2389150, 0.2510757, 0.9108029, 1.0566703, 0.2017566, 0.4542836, 1.0635974, 0.3674849, 0.1636298 | 7.762017, 7.537917, 7.451294, 8.303507, 5.789990, 6.708877, 6.950702, 7.079058, 8.038645, 7.425258, 7.743959, 7.729417, 6.226316, 6.221490, 6.068456, 7.161888, 7.567652, 6.010780, 7.687761, 6.452530 | 0.00e+00 | welch_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| NB | CCDC64 | 4.250962, 5.667639, 5.997046, 5.333400, 5.722461, 4.979527, 5.579268, 6.195898, 5.974999, 5.705686, 5.910760, 4.782449, 5.346567, 4.908316, 6.692405 | 3.228049, 3.698218, 3.001802, 3.533563, 3.471187, 3.195348, 3.549669, 2.400538, 1.669027, 2.643856, 4.671859, 3.381283, 3.638074, 3.125982, 1.773996 | 0.00e+00 | student_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| NB | HIST1H1B | 0.04278404, 0.73994064, 0.13763218, 0.29878223, 0.51612818, 0.11116746, 0.55591464, 0.00000000, 0.45428360, 0.66001219, 0.18916444, 0.07052473, 0.21424465, 0.38968015, 0.09774809 | 4.866413, 4.683696, 5.180307, 6.172528, 4.809929, 3.408712, 3.589763, 7.021924, 5.826040, 3.648465, 4.938756, 6.137913, 3.875780, 2.378512, 3.153805 | 0.00e+00 | welch_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| SS | HIST1H1B | 0.18903382, 0.05672351, 0.48542683, 0.07038933, 0.00000000, 0.21424463, 0.95605665, 0.08406426, 0.08420307, 0.42223300 | 5.480911, 4.917432, 5.080231, 4.204767, 3.982765, 5.044831, 5.594549, 3.760221, 5.142413, 1.469886, 6.399171, 6.007644, 4.563158, 6.052677, 2.558268 | 3.59e-05 | mann_whitney | \*\*\*\* | 3.59e-05 | \*\*\*\* |
| SS | MTMR2 | 4.370164, 4.472532, 4.458776, 4.869871, 4.771357, 4.832883, 4.823749, 5.116448, 5.057450, 4.797532, 4.689858, 4.841973, 3.653083, 4.622930, 4.927896 | 3.344828, 3.496974, 3.235727, 3.217231, 3.442280, 2.914565, 2.435629, 2.639232, 2.992768, 2.965323, 1.903038, 2.811471, 2.954196, 2.757023, 1.220330 | 3.40e-06 | mann_whitney | \*\*\*\* | 6.80e-06 | \*\*\*\* |
| WT | HIST1H1B | 1.7355312, 0.3449496, 1.4751447, 1.5110093, 1.3449054, 1.5310872, 0.5559146, 1.7398738, 1.9411447, 2.0976171, 1.5509181, 2.0250578, 1.6508209, 1.7866591 | 7.276962, 7.303141, 8.023810, 7.387845, 7.817879, 7.113117, 8.372647, 7.557042, 8.135453, 7.485185, 8.715585, 7.308248, 6.970509, 7.379292, 8.057342 | 5.10e-06 | mann_whitney | \*\*\*\* | 5.10e-06 | \*\*\*\* |
| WT | SYPL1 | 5.905253, 5.808943, 6.436934, 6.170353, 5.358566, 6.251149, 6.023431, 5.721971, 7.211061, 5.536103, 6.359660, 5.872318, 6.028748, 6.511882, 6.706369 | 4.025029, 3.099295, 3.587365, 3.084064, 2.935460, 3.132577, 3.169925, 3.716991, 2.669027, 4.016140, 3.092546, 3.755956, 4.202418, 3.788686, 4.265287 | 0.00e+00 | student_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |
| aRMS | HIST1H1B | 0.13750352, 0.65076456, 0.59454855, 2.45680615, 0.23878686, 0.09761080, 0.08406426, 1.46988598, 0.94110631, 0.54596837, 1.01435529, 0.41142625, 0.16349873, 0.22650853 | 6.483655, 8.466219, 5.283551, 9.088550, 6.588265, 7.465158, 7.260120, 6.940519, 4.895303, 8.416375, 6.124121, 4.242603, 7.101398, 5.907852, 8.798115 | 5.10e-06 | mann_whitney | \*\*\*\* | 5.10e-06 | \*\*\*\* |
| aRMS | TM2D1 | 5.010780, 4.244126, 4.599318, 4.748461, 5.377124, 5.503985, 5.365623, 5.294988, 4.625270, 5.261531, 4.914565, 4.349082, 5.029011, 4.154616, 5.081936 | 3.590961, 2.538538, 3.340562, 3.482848, 3.928844, 3.292782, 2.691534, 3.704872, 3.065228, 4.569856, 2.931683, 2.963474, 2.531069, 2.957915, 2.965323 | 0.00e+00 | student_t | \*\*\*\* | 0.00e+00 | \*\*\*\* |

Apply stars to figures

``` r
merged_plots <- map(diseases, function(d) {

  median_gene <- median_gene_map[[d]]

  # HIST1H1B expression
  df_h1b <- combined_long %>%
    filter(Disease == d, Gene == "HIST1H1B") %>%
    mutate(PlotType = "HIST1H1B (PolyA-)")

  # Median-ratio gene expression
  df_median <- combined_long %>%
    filter(Disease == d, Gene == median_gene) %>%
    mutate(PlotType = "Median Ratio Gene (PolyA+)")

  # Combine
  df_combined <- bind_rows(df_h1b, df_median)

  # SIGNIFICANCE STARS
  star_df <- tibble(
    Disease = d,
    Gene = c("HIST1H1B", median_gene),
    PlotType = c("HIST1H1B (PolyA-)", "Median Ratio Gene (PolyA+)"),
    stars_adj = sig_results %>%
      filter(Disease == d, Gene %in% c("HIST1H1B", median_gene)) %>%
      arrange(match(Gene, c("HIST1H1B", median_gene))) %>%
      pull(stars_adj)
  )

  # y-position for stars (slightly above max)
  y_star <- max(df_combined$Expression, na.rm = TRUE) * 1.05

  ggplot(df_combined, aes(x = Gene, y = Expression, fill = Compendia)) +
    geom_boxplot(
      alpha = 0.7,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      linewidth = 0.4
    ) +
    geom_jitter(
      aes(color = Compendia),
      position = position_jitterdodge(jitter.width = 0.2,
                                      dodge.width = 0.8),
      size = 1.2,
      alpha = 0.6
    ) +
    # ADD STARS
    geom_text(
      data = star_df,
      aes(x = Gene, y = y_star, label = stars_adj),
      inherit.aes = FALSE,
      size = 6
    ) +
    scale_fill_compendia() +
    scale_color_compendia() +
    coord_cartesian(ylim = c(0, y_star * 1.1)) +
    facet_wrap(~ PlotType, ncol = 2, scales = "free_x") +
    labs(
      title = paste("Expression of Representative PolyA+ and PolyA- Gene in", d),
      x = "Gene",
      y = "Expression log2(TPM+1)"
    ) +
    color_theme()
})
merged_plots
```

    [[1]]

![](Fig1E_files/figure-commonmark/Fig1E-1.png)


    [[2]]

![](Fig1E_files/figure-commonmark/Fig1E-2.png)


    [[3]]

![](Fig1E_files/figure-commonmark/Fig1E-3.png)


    [[4]]

![](Fig1E_files/figure-commonmark/Fig1E-4.png)


    [[5]]

![](Fig1E_files/figure-commonmark/Fig1E-5.png)


    [[6]]

![](Fig1E_files/figure-commonmark/Fig1E-6.png)

## Fig 1F

Histogram of median ratio distributions, calculated as the median
expression level in PolyA samples divided by the median expression level
in RiboD samples

``` r
ss_ratios <- SS_gene_medians %>%
  select(Gene, Ratio = SS_median_ratio) %>%
  mutate(Disease = "SS")

arms_ratios <- aRMS_gene_medians %>%
  select(Gene = Gene, Ratio = aRMS_median_ratio) %>%
  mutate(Disease = "aRMS")

wt_ratios <- WT_gene_medians %>%
  select(Gene, Ratio = WT_median_ratio) %>%
  mutate(Disease = "WT")

nb_ratios <- NB_gene_medians %>%
  select(Gene = Gene, Ratio = NB_median_ratio) %>%
  mutate(Disease = "NB")

all_ratios <- ALL_gene_medians %>%
  select(Gene = Gene, Ratio = ALL_median_ratio) %>%
  mutate(Disease = "ALL")

aml_ratios <- AML_gene_medians %>%
  select(Gene = Gene, Ratio = AML_median_ratio) %>%
  mutate(Disease = "AML")

#Combining median ratios
ratio_df <- bind_rows(ss_ratios, arms_ratios, wt_ratios, nb_ratios, aml_ratios, all_ratios)
```

``` r
theme_Fig1F <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      plot.title = element_text(face = "bold", size = base_size * 1.1, hjust = 0.5),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 15), 
      plot.margin = margin(10, 10, 10, 10)
    )
}
```

``` r
SS_HIST1H1B <- SS_gene_medians_hugo %>%
  filter(Gene %in% c("HIST1H1B"))

aRMS_HIST1H1B <- aRMS_gene_median_hugo %>%
  filter(Gene %in% c("HIST1H1B"))

WT_HIST1H1B <- WT_gene_medians_hugo %>%
  filter(Gene %in% c("HIST1H1B"))

NB_HIST1H1B <- NB_gene_medians_hugo %>%
  filter(Gene %in% c("HIST1H1B"))

AML_HIST1H1B <- AML_gene_medians_hugo %>%
  filter(Gene %in% c("HIST1H1B"))

ALL_HIST1H1B <- ALL_gene_medians_hugo %>%
  filter(Gene %in% c("HIST1H1B"))
```

``` r
Fig1F_aRMS <- ggplot(arms_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = aRMS_most_average_gene_hugo$aRMS_median_ratio, color = "#0072B2", linetype = "dashed") +
  geom_text(aes(x = 2.2, y = 2000, label = aRMS_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  geom_vline(xintercept = aRMS_HIST1H1B$aRMS_median_ratio, color = "#E69F00", linetype = "dashed") +
  geom_text(aes(x = 0.8, y = 3000, label = aRMS_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in aRMS Samples",
    x = "aRMS Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_aRMS
```

    Warning in geom_text(aes(x = 2.2, y = 2000, label = aRMS_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.8, y = 3000, label = aRMS_HIST1H1B$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 27707 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig1F%20aRMS-1.png)

``` r
Fig1F_SS <- ggplot(ss_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = SS_most_average_gene_hugo$SS_median_ratio, color = "#0072B2", linetype = "dashed") +
  
  geom_text(aes(x = 2.2, y = 2000, label = SS_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  
  geom_vline(xintercept = SS_HIST1H1B$SS_median_ratio, color = "#E69F00", linetype = "dashed") +
  
  geom_text(aes(x = 0.8, y = 3000, label = SS_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in SS Samples",
    x = "SS Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_SS
```

    Warning in geom_text(aes(x = 2.2, y = 2000, label = SS_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.8, y = 3000, label = SS_HIST1H1B$Gene), angle = 0, : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 29257 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig1F%20SS-1.png)

``` r
Fig1F_WT <- ggplot(wt_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = WT_most_average_gene_hugo$WT_median_ratio, color = "#0072B2", linetype = "dashed") +
  
  geom_text(aes(x = 2.2, y = 2000, label = WT_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  
  geom_vline(xintercept = WT_HIST1H1B$WT_median_ratio, color = "#E69F00", linetype = "dashed") +
  
  geom_text(aes(x = 0.95, y = 3000, label = WT_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in WT Samples",
    x = "WT Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_WT
```

    Warning in geom_text(aes(x = 2.2, y = 2000, label = WT_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.95, y = 3000, label = WT_HIST1H1B$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 26797 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig%201f%20WT-1.png)

``` r
Fig1F_NB <- ggplot(nb_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = NB_most_average_gene_hugo$NB_median_ratio, color = "#0072B2", linetype = "dashed") +
  
  geom_text(aes(x = 2.4, y = 2000, label = NB_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  
  geom_vline(xintercept = NB_HIST1H1B$NB_median_ratio, color = "#E69F00", linetype = "dashed") +
  
  geom_text(aes(x = 0.9, y = 3000, label = NB_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in NB Samples",
    x = "NB Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_NB
```

    Warning in geom_text(aes(x = 2.4, y = 2000, label = NB_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.9, y = 3000, label = NB_HIST1H1B$Gene), angle = 0, : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 25112 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig1F%20NB-1.png)

``` r
Fig1F_ALL <- ggplot(all_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = ALL_most_average_gene_hugo$ALL_median_ratio, color = "#0072B2", linetype = "dashed") +
  
  geom_text(aes(x = 2.9, y = 2000, label = ALL_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  
  geom_vline(xintercept = ALL_HIST1H1B$ALL_median_ratio, color = "#E69F00", linetype = "dashed") +
  
  geom_text(aes(x = 0.9, y = 2500, label = ALL_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in ALL Samples",
    x = "ALL Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_ALL
```

    Warning in geom_text(aes(x = 2.9, y = 2000, label = ALL_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.9, y = 2500, label = ALL_HIST1H1B$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 31150 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig1F%20ALL-1.png)

``` r
Fig1F_AML <- ggplot(aml_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = AML_most_average_gene_hugo$AML_median_ratio, color = "#0072B2", linetype = "dashed") +
  
  geom_text(aes(x = 2.5, y = 1500, label = AML_most_average_gene_hugo$Gene), angle = 0, vjust = -0.5, color = "#0072B2", size = 4) +
  
  geom_vline(xintercept = AML_HIST1H1B$AML_median_ratio, color = "#E69F00", linetype = "dashed") +
  
  geom_text(aes(x = 0.9, y = 1500, label = AML_HIST1H1B$Gene), angle = 0, vjust = -0.5, color = "#E69F00", size = 4) +
  
  coord_cartesian(xlim = c(0, 10)) +
  labs(
#    title = "Distribution of polyA / riboD Median Ratios in AML Samples",
    x = "AML Median Expression Ratio (polyA / riboD)",
    y = "Number of Genes"
  ) +
  theme_Fig1F()

Fig1F_AML
```

    Warning in geom_text(aes(x = 2.5, y = 1500, label = AML_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.9, y = 1500, label = AML_HIST1H1B$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 30953 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1E_files/figure-commonmark/fig1F%20AML-1.png)

Session Info

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
     date     2026-06-12
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     ! package      * version date (UTC) lib source
     P bit            4.6.0   2025-03-06 [?] RSPM
     P bit64          4.8.0   2026-04-21 [?] RSPM
     P cli            3.6.5   2025-04-23 [?] RSPM
     P crayon         1.5.3   2024-06-20 [?] RSPM
     P digest         0.6.37  2024-08-19 [?] RSPM
     P dplyr        * 1.2.1   2026-04-03 [?] RSPM
     P evaluate       1.0.5   2025-08-27 [?] RSPM
     P farver         2.1.2   2024-05-13 [?] RSPM
     P fastmap        1.2.0   2024-05-15 [?] RSPM
     P forcats      * 1.0.1   2025-09-25 [?] RSPM
     P generics       0.1.4   2025-05-09 [?] RSPM
     P ggplot2      * 4.0.3   2026-04-22 [?] RSPM
     P glue           1.8.0   2024-09-30 [?] RSPM
     P gtable         0.3.6   2024-10-25 [?] RSPM
     P hms            1.1.4   2025-10-17 [?] RSPM
     P htmltools      0.5.8.1 2024-04-04 [?] RSPM
     P jsonlite       2.0.0   2025-03-27 [?] RSPM
     P knitr          1.50    2025-03-16 [?] RSPM
     P labeling       0.4.3   2023-08-29 [?] RSPM
     P lifecycle      1.0.5   2026-01-08 [?] RSPM
     P lubridate    * 1.9.5   2026-02-04 [?] RSPM
     P magrittr       2.0.5   2026-04-04 [?] RSPM
     P pillar         1.11.1  2025-09-17 [?] RSPM
     P pkgconfig      2.0.3   2019-09-22 [?] RSPM
     P purrr        * 1.2.2   2026-04-10 [?] RSPM
     P R6             2.6.1   2025-02-15 [?] RSPM
     P RColorBrewer   1.1-3   2022-04-03 [?] RSPM
     P readr        * 2.2.0   2026-02-19 [?] RSPM
     P rlang          1.2.0   2026-04-06 [?] RSPM
     P rmarkdown      2.30    2025-09-28 [?] RSPM
     P rstudioapi     0.18.0  2026-01-16 [?] RSPM
     P S7             0.2.2   2026-04-22 [?] RSPM
     P scales         1.4.0   2025-04-24 [?] RSPM
     P sessioninfo    1.2.3   2025-02-05 [?] CRAN (R 4.5.0)
     P stringi        1.8.7   2025-03-27 [?] RSPM
     P stringr      * 1.6.0   2025-11-04 [?] RSPM
     P tibble       * 3.3.1   2026-01-11 [?] RSPM
     P tidyr        * 1.3.2   2025-12-19 [?] RSPM
     P tidyselect     1.2.1   2024-03-11 [?] RSPM
     P tidyverse    * 2.0.0   2023-02-22 [?] RSPM
     P timechange     0.4.0   2026-01-29 [?] RSPM
     P tzdb           0.5.0   2025-03-15 [?] RSPM
     P vctrs          0.7.3   2026-04-11 [?] RSPM
     P vroom          1.7.1   2026-03-31 [?] RSPM
     P withr          3.0.2   2024-10-28 [?] RSPM
     P xfun           0.55    2025-12-16 [?] CRAN (R 4.5.2)
     P yaml           2.3.10  2024-07-26 [?] RSPM

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/Fig_1/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
