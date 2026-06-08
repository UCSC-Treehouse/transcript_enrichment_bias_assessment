# correlation_analysis_TEBA
Maryke Grobler

# correlation_analysis

Maryke Grobler

## Correlation analysis

The goal of this analysis is to determine which samples are the most
highly correlated within each disease type, then subset the gene
expression to only include those samples. The subsetted expression file
will be used for downstream differential expression analysis.

This is using gene expression from the Treehouse Tumor v25 RiboD and
PolyA compendia. Diseases include synovial sarcoma (SS), alveolar
rhabdomyosarcoma (aRMS), Wilms tumor (WT), neuroblastoma (NB), acute
myeloid leukemia (AML), and acute lymphoblastic leukemia (ALL).

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
library(corrr)
```

``` r
polyA_clinical <- read_tsv("https://xena.treehouse.gi.ucsc.edu:443/download/clinical_Treehouse-Tumor-Compendium-25.01-PolyA_20250131v1.tsv")
```

    Rows: 13359 Columns: 12
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (12): th_dataset_id, age_at_dx, disease, icd_disease, organism, pedaya, ...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
riboD_clinical <- read_tsv("https://xena.treehouse.gi.ucsc.edu:443/download/clinical_Treehouse-Tumor-Compendium-25.01-RiboD_20250306v1.tsv")
```

    Rows: 2079 Columns: 12
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (12): th_dataset_id, age_at_dx, disease, icd_disease, organism, pedaya, ...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# create column in each file that specifies with compendia the file is from
riboD_clinical_lib <- riboD_clinical %>%
  mutate(lib_prep="riboD")

polyA_clinical_lib <- polyA_clinical %>%
  mutate(lib_prep="polyA")

# merge clinical files
clin_merged <- bind_rows(polyA_clinical_lib, riboD_clinical_lib)
```

## SS and aRMS

``` r
expression_SS_aRMS <- read_tsv("../input_data/expression_SS_aRMS.tsv.gz", show_col_types = FALSE)
```

``` r
SS_aRMS <- c("synovial sarcoma", "alveolar rhabdomyosarcoma")

clin_merged_SS_aRMS <- clin_merged %>%
  filter(disease %in% c(SS_aRMS)) %>%
  group_by(disease) %>%
  mutate(disease = case_when(
    disease == "synovial sarcoma" ~ "SS",
    disease == "alveolar rhabdomyosarcoma" ~ "aRMS"
    )) %>%
  unite(col = "disease_and_prep", disease, lib_prep, sep = "_") %>%
  select(th_dataset_id, disease_and_prep)

#write_tsv(clin_merged_SS_aRMS, "../data/clin_merged_SS_aRMS.tsv")
```

The datasets include some samples from the same patient, which may be
skew the correlation scores since those samples would be more highly
correlated to each other than any of the other samples.

``` r
# filtering out samples that contain S02 or S03
# in the clinical file
clin_SS_aRMS_dups_removed <- clin_merged_SS_aRMS %>%
  filter(!grepl("S02$", th_dataset_id),
         !grepl("S03$", th_dataset_id))

# in the expression file
expr_SS_aRMS_dups_removed <- expression_SS_aRMS %>%
  pivot_longer(cols = -Gene,
               names_to = "THID",
               values_to = "est_count") %>%
  filter(!grepl("S02$", THID),
         !grepl("S03$", THID)) %>%
  pivot_wider(names_from = "THID",
              values_from = "est_count")
```

``` r
# the first column in expr_SS_aRMS_dups_removed is "Gene", but the gene names should be the names of each row, not as a separate column.

cleaned_expr_SS_aRMS_dups_rem <- expr_SS_aRMS_dups_removed %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")
```

``` r
# filtering out genes with a count less than 10 across all samples
filt_expr_SS_aRMS_dups_rem <- cleaned_expr_SS_aRMS_dups_rem %>%
  filter(if_all(everything(), ~.x > 10))
```

``` r
# filtering clinical file 
# DR means Duplicate samples Removed
clin_SS_polyA_DR <- clin_SS_aRMS_dups_removed %>%
  filter(disease_and_prep == "SS_polyA")

clin_SS_riboD_DR <- clin_SS_aRMS_dups_removed %>%
  filter(disease_and_prep == "SS_riboD")

clin_aRMS_polyA_DR <- clin_SS_aRMS_dups_removed %>%
  filter(disease_and_prep == "aRMS_polyA")

clin_aRMS_riboD_DR <- clin_SS_aRMS_dups_removed %>%
  filter(disease_and_prep == "aRMS_riboD")

# filtering expression file
expression_SS_polyA_DR <- filt_expr_SS_aRMS_dups_rem %>%
  select(all_of(c(clin_SS_polyA_DR$th_dataset_id)))

expression_SS_riboD_DR <- filt_expr_SS_aRMS_dups_rem %>%
  select(all_of(c(clin_SS_riboD_DR$th_dataset_id)))

expression_aRMS_polyA_DR <- filt_expr_SS_aRMS_dups_rem %>%
  select(all_of(c(clin_aRMS_polyA_DR$th_dataset_id)))

expression_aRMS_riboD_DR <- filt_expr_SS_aRMS_dups_rem %>%
  select(all_of(c(clin_aRMS_riboD_DR$th_dataset_id)))
```

``` r
cor_SS_polyA_DR <- correlate(expression_SS_polyA_DR, method = "spearman", quiet = TRUE)
cor_SS_riboD_DR <- correlate(expression_SS_riboD_DR, method = "spearman", quiet = TRUE)
cor_aRMS_polyA_DR <- correlate(expression_aRMS_polyA_DR, method = "spearman", quiet = TRUE)
cor_aRMS_riboD_DR <- correlate(expression_aRMS_riboD_DR, method = "spearman", quiet = TRUE)
```

``` r
# SS_polyA
sum_cor_SS_polyA_DR <- cor_SS_polyA_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_SS_polyA_DR, n = 5)
```

| term            | min_value | max_value | median_value | mean_value |     stdev |
|:----------------|----------:|----------:|-------------:|-----------:|----------:|
| THR39_1373_S01  | 0.5761078 | 0.9294128 |    0.8860025 |  0.8551305 | 0.0754844 |
| TCGA-WK-A8XT-01 | 0.6045318 | 0.9239853 |    0.8784215 |  0.8534203 | 0.0702340 |
| TH40_2281_S01   | 0.5468642 | 0.9283242 |    0.8781524 |  0.8477257 | 0.0797938 |
| THR39_1375_S01  | 0.5277929 | 0.9719540 |    0.8780302 |  0.8490935 | 0.0916370 |
| THR31_1392_S01  | 0.5368382 | 0.9329047 |    0.8740251 |  0.8403389 | 0.0906675 |

``` r
# SS_riboD
sum_cor_SS_riboD_DR <- cor_SS_riboD_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise( #calculate the stats
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_SS_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR51_4556_S01 | 0.6205281 | 0.9383515 |    0.8937483 |  0.8542862 | 0.0832300 |
| THR51_4558_S01 | 0.5964468 | 0.9491028 |    0.8837590 |  0.8486666 | 0.0943203 |
| THR51_4554_S01 | 0.5979186 | 0.9383515 |    0.8795101 |  0.8397882 | 0.0875092 |
| THR24_3992_S01 | 0.6337958 | 0.9446409 |    0.8707665 |  0.8414040 | 0.0839337 |
| THR24_3887_S01 | 0.6067254 | 0.9491028 |    0.8625286 |  0.8361579 | 0.0930422 |

``` r
# aRMS_polyA
# convert to long format, reshapes the data from wide (74 columns) to long format (one row per sample pair)
sum_cor_aRMS_polyA_DR <- cor_aRMS_polyA_DR %>%
  pivot_longer(cols = -term,
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_aRMS_polyA_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR29_0788_S01 | 0.6672448 | 0.9381194 |    0.8692229 |  0.8615073 | 0.0483446 |
| THR29_0775_S01 | 0.6475805 | 0.9290565 |    0.8680855 |  0.8607676 | 0.0488149 |
| THR29_0757_S01 | 0.6663602 | 0.9315209 |    0.8621625 |  0.8540454 | 0.0514703 |
| THR29_0762_S01 | 0.6756412 | 0.9347398 |    0.8598925 |  0.8506785 | 0.0482096 |
| THR29_0771_S01 | 0.6311285 | 0.9286304 |    0.8588158 |  0.8512309 | 0.0521252 |

``` r
# aRMS_riboD
# convert to long format, reshapes the data from wide (74 columns) to long format (one row per sample pair)
sum_cor_aRMS_riboD_DR <- cor_aRMS_riboD_DR %>%
  pivot_longer(cols = -term,
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_aRMS_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_3244_S01 | 0.6200743 | 0.9358870 |    0.8835668 |  0.8507649 | 0.0904862 |
| THR24_3371_S01 | 0.6046084 | 0.9367303 |    0.8793823 |  0.8399199 | 0.0867619 |
| THR24_3181_S01 | 0.6323477 | 0.9464215 |    0.8776794 |  0.8509646 | 0.0875724 |
| THR24_3178_S01 | 0.5859419 | 0.9367303 |    0.8771115 |  0.8380068 | 0.0886288 |
| THR24_3606_S01 | 0.6057197 | 0.9369294 |    0.8750470 |  0.8395270 | 0.0975532 |

``` r
# this was a sample that was only 18th out of 20 most correlated samples, but it clustered together with the other SS riboD samples, so I'm including it
THR49_SS_riboD <- sum_cor_SS_riboD_DR %>%
  filter(term == "THR49_4787_S01") 

filtered_SS_riboD_list <- sum_cor_SS_riboD_DR %>%
  filter(term != "TH38_1433_S01") %>% # this was 10 out of 20 in the most correlated, but it clustered in the SS polyA samples, so I'm removing it
  slice_head(n = 14) %>% #taking the top most correlated samples
  rbind(THR49_SS_riboD) %>% #adding this THID 
  select(term) %>%
  mutate(disease_and_prep = "SS_riboD") # need to add back in the category

filtered_SS_polyA_list <- sum_cor_SS_polyA_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "SS_polyA")

filtered_aRMS_polyA_list <- sum_cor_aRMS_polyA_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "aRMS_polyA")

filtered_aRMS_riboD_list <- sum_cor_aRMS_riboD_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "aRMS_riboD")
```

``` r
filtered_SS_aRMS_list <- bind_rows(filtered_SS_riboD_list, filtered_SS_polyA_list, filtered_aRMS_riboD_list, filtered_aRMS_polyA_list)

write_tsv(filtered_SS_aRMS_list, "../input_data/filtered_SS_aRMS_list.tsv")
```

## WT and NB

``` r
expression_WT_NB <- read_tsv("../input_data/expression_WT_NB.tsv.gz", show_col_types = FALSE)
```

``` r
WT_NB <- c("wilms tumor", "neuroblastoma")

clin_merged_WT_NB <- clin_merged %>%
  filter(disease %in% c(WT_NB)) %>%
  group_by(disease) %>%
  mutate(disease = case_when(
    disease == "wilms tumor" ~ "WT",
    disease == "neuroblastoma" ~ "NB"
    )) %>%
  unite(col = "disease_and_prep", disease, lib_prep, sep = "_") %>%
  select(th_dataset_id, disease_and_prep) 

#write_tsv(clin_merged_WT_NB, "../data/clin_merged_WT_NB.tsv")
```

``` r
# filtering out samples that contain S02, S03, S04
# in the clinical file
clin_WT_NB_dups_removed <- clin_merged_WT_NB %>%
  filter(!grepl("S02$", th_dataset_id), # the $ in a regular expression indicates the end of a line
         !grepl("S03$", th_dataset_id),
         !grepl("S04$", th_dataset_id))

# in the expression file
expr_WT_NB_dups_removed <- expression_WT_NB %>%
  pivot_longer(cols = -Gene,
               names_to = "THID",
               values_to = "est_count") %>%
  filter(!grepl("S02$", THID),
         !grepl("S03$", THID),
         !grepl("S04$", THID)) %>%
  pivot_wider(names_from = "THID",
              values_from = "est_count")
```

``` r
# the first column in expr_SS_aRMS_dups_removed is "Gene", but the gene names should be the names of each row, not as a separate column.

cleaned_expr_WT_NB_dups_removed <- expr_WT_NB_dups_removed %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")
```

``` r
filt_expr_WT_NB_dups_rem <- cleaned_expr_WT_NB_dups_removed %>%
  filter(if_all(everything(), ~.x > 10))
```

``` r
# filtering clinical file 
# DR means Duplicate samples Removed
clin_WT_polyA_DR <- clin_WT_NB_dups_removed %>%
  filter(disease_and_prep == "WT_polyA")

clin_WT_riboD_DR <- clin_WT_NB_dups_removed %>%
  filter(disease_and_prep == "WT_riboD")

clin_NB_polyA_DR <- clin_WT_NB_dups_removed %>%
  filter(disease_and_prep == "NB_polyA")

clin_NB_riboD_DR <- clin_WT_NB_dups_removed %>%
  filter(disease_and_prep == "NB_riboD")

# filtering expression file
expression_WT_polyA_DR <- filt_expr_WT_NB_dups_rem %>%
  select(all_of(c(clin_WT_polyA_DR$th_dataset_id)))

expression_WT_riboD_DR <- filt_expr_WT_NB_dups_rem %>%
  select(all_of(c(clin_WT_riboD_DR$th_dataset_id)))

expression_NB_polyA_DR <- filt_expr_WT_NB_dups_rem %>%
  select(all_of(c(clin_NB_polyA_DR$th_dataset_id)))

expression_NB_riboD_DR <- filt_expr_WT_NB_dups_rem %>%
  select(all_of(c(clin_NB_riboD_DR$th_dataset_id)))
```

``` r
# Run correlate()
cor_WT_polyA_DR <- correlate(expression_WT_polyA_DR, method = "spearman", quiet = TRUE)
cor_WT_riboD_DR <- correlate(expression_WT_riboD_DR, method = "spearman", quiet = TRUE)
cor_NB_polyA_DR <- correlate(expression_NB_polyA_DR, method = "spearman", quiet = TRUE)
cor_NB_riboD_DR <- correlate(expression_NB_riboD_DR, method = "spearman", quiet = TRUE)
```

``` r
# WT_polyA
sum_cor_WT_polyA_DR <- cor_WT_polyA_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_WT_polyA_DR, n = 5)
```

| term                | min_value | max_value | median_value | mean_value |     stdev |
|:--------------------|----------:|----------:|-------------:|-----------:|----------:|
| TARGET-50-PAJNCZ-01 | 0.5998466 | 0.9695092 |    0.9152583 |  0.8945627 | 0.0695764 |
| TARGET-50-PAEBXA-01 | 0.6667870 | 0.9559104 |    0.9135948 |  0.8965815 | 0.0547466 |
| TARGET-50-PALERC-01 | 0.6021529 | 0.9765217 |    0.9128599 |  0.8956669 | 0.0630733 |
| TARGET-50-PALGAZ-01 | 0.6166671 | 0.9566575 |    0.9128436 |  0.8916217 | 0.0653862 |
| TARGET-50-PAJLNJ-01 | 0.6327158 | 0.9710931 |    0.9118714 |  0.8940442 | 0.0585125 |

``` r
# WT_riboD
sum_cor_WT_riboD_DR <- cor_WT_riboD_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise( #calculate the stats
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_WT_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_3218_S01 | 0.7465149 | 0.9440634 |    0.8899325 |  0.8708309 | 0.0558929 |
| THR24_4194_S01 | 0.6699285 | 0.9398217 |    0.8820968 |  0.8468491 | 0.0814571 |
| THR24_4284_S01 | 0.7260917 | 0.9393614 |    0.8820257 |  0.8563592 | 0.0664292 |
| THR24_4369_S01 | 0.6892870 | 0.9393614 |    0.8819004 |  0.8435356 | 0.0790015 |
| THR24_4226_S01 | 0.7217054 | 0.9440634 |    0.8811747 |  0.8512656 | 0.0688677 |

``` r
# NB_polyA
sum_cor_NB_polyA_DR <- cor_NB_polyA_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_NB_polyA_DR, n = 5)
```

| term                | min_value | max_value | median_value | mean_value |     stdev |
|:--------------------|----------:|----------:|-------------:|-----------:|----------:|
| TARGET-30-PASUML-01 | 0.5723757 | 0.9482206 |    0.8615253 |  0.8457083 | 0.0640410 |
| TARGET-30-PASEGA-01 | 0.5596472 | 0.9109113 |    0.8548460 |  0.8334860 | 0.0637989 |
| TARGET-30-PAPUAR-01 | 0.6106996 | 0.9566123 |    0.8542641 |  0.8400699 | 0.0707671 |
| TARGET-30-PATHVK-01 | 0.5670192 | 0.9170777 |    0.8534949 |  0.8303141 | 0.0731170 |
| TARGET-30-PATFCY-01 | 0.5601384 | 0.9354618 |    0.8522847 |  0.8264972 | 0.0759165 |

``` r
# NB_riboD
sum_cor_NB_riboD_DR <- cor_NB_riboD_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise( #calculate the stats
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_NB_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_4310_S01 | 0.5683897 | 0.9506633 |    0.8739380 |  0.8518428 | 0.0768598 |
| THR24_2779_S01 | 0.5719067 | 0.9312754 |    0.8656171 |  0.8382691 | 0.0730688 |
| THR24_3516_S01 | 0.5290457 | 0.9506633 |    0.8617932 |  0.8390842 | 0.0869628 |
| THR24_4114_S01 | 0.5566618 | 0.9447870 |    0.8598830 |  0.8369181 | 0.0850224 |
| THR24_4137_S01 | 0.5454747 | 0.9331739 |    0.8598216 |  0.8348221 | 0.0779318 |

The top 15 most correlated samples were all in the same cluster in the
tumormap, so there was minimal issue in creating the subsetted sample
IDs.

``` r
filtered_WT_polyA_list <- sum_cor_WT_polyA_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "WT_polyA")

filtered_WT_riboD_list <- sum_cor_WT_riboD_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "WT_riboD")

filtered_NB_polyA_list <- sum_cor_NB_polyA_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "NB_polyA")

filtered_NB_riboD_list <- sum_cor_NB_riboD_DR %>%
  slice_head(n = 15) %>%
  select(term) %>%
  mutate(disease_and_prep = "NB_riboD")
```

``` r
filtered_WT_NB_list <- bind_rows(filtered_WT_riboD_list, filtered_WT_polyA_list, filtered_NB_riboD_list, filtered_NB_polyA_list)

write_tsv(filtered_WT_NB_list, "../input_data/filtered_WT_NB_list.tsv")
```

## ALL and AML

``` r
expression_ALL_AML <- read_tsv("../input_data/expression_ALL_AML.tsv.gz", show_col_types = FALSE)
```

``` r
ALL_AML <- c("acute lymphoblastic leukemia", "acute myeloid leukemia")

clin_merged_ALL_AML <- clin_merged %>%
  filter(disease %in% c(ALL_AML)) %>%
  group_by(disease) %>%
  mutate(disease = case_when(
    disease == "acute myeloid leukemia" ~ "AML",
    disease == "acute lymphoblastic leukemia" ~ "ALL"
    )) %>%
  unite(col = "disease_and_prep", disease, lib_prep, sep = "_") %>%
  select(th_dataset_id, disease_and_prep)

#write_tsv(clin_merged_ALL_AML, "../data/clin_merged_ALL_AML.tsv")
```

``` r
# in the clinical file
clin_ALL_AML_noTARGET_dups_removed <- clin_merged_ALL_AML %>%
  filter(!grepl("TARGET", th_dataset_id),
         !grepl("S02$", th_dataset_id),
         !grepl("S03$", th_dataset_id),
         !grepl("S04$", th_dataset_id))

# in the expression file
expr_ALL_AML_noTARGET_dups_removed <- expression_ALL_AML %>%
  pivot_longer(cols = -Gene,
               names_to = "THID",
               values_to = "est_count") %>%
  filter(!grepl("TARGET", THID),
         !grepl("S02$", THID),
         !grepl("S03$", THID),
         !grepl("S04$", THID)) %>%
  pivot_wider(names_from = "THID",
              values_from = "est_count")
```

``` r
cleaned_expr_ALL_AML_noTARGET_DR <- expr_ALL_AML_noTARGET_dups_removed %>%
  remove_rownames %>%
  column_to_rownames(var = "Gene")
```

``` r
filt_expr_ALL_AML_dups_rem <- cleaned_expr_ALL_AML_noTARGET_DR %>%
  filter(if_all(everything(), ~.x > 10))
```

``` r
# filtering clinical file 
# DR means Duplicate samples Removed
clin_ALL_polyA_DR <- clin_ALL_AML_noTARGET_dups_removed %>%
  filter(disease_and_prep == "ALL_polyA")

clin_ALL_riboD_DR <- clin_ALL_AML_noTARGET_dups_removed %>%
  filter(disease_and_prep == "ALL_riboD")

clin_AML_polyA_DR <- clin_ALL_AML_noTARGET_dups_removed %>%
  filter(disease_and_prep == "AML_polyA")

clin_AML_riboD_DR <- clin_ALL_AML_noTARGET_dups_removed %>%
  filter(disease_and_prep == "AML_riboD")

# filtering expression file
expression_ALL_polyA_DR <- filt_expr_ALL_AML_dups_rem %>%
  select(all_of(c(clin_ALL_polyA_DR$th_dataset_id)))

expression_ALL_riboD_DR <- filt_expr_ALL_AML_dups_rem %>%
  select(all_of(c(clin_ALL_riboD_DR$th_dataset_id)))

expression_AML_polyA_DR <- filt_expr_ALL_AML_dups_rem %>%
  select(all_of(c(clin_AML_polyA_DR$th_dataset_id)))

expression_AML_riboD_DR <- filt_expr_ALL_AML_dups_rem %>%
  select(all_of(c(clin_AML_riboD_DR$th_dataset_id)))
```

``` r
cor_ALL_polyA_DR <- correlate(expression_ALL_polyA_DR, method = "spearman", quiet = TRUE)
cor_ALL_riboD_DR <- correlate(expression_ALL_riboD_DR, method = "spearman", quiet = TRUE)
cor_AML_polyA_DR <- correlate(expression_AML_polyA_DR, method = "spearman", quiet = TRUE)
cor_AML_riboD_DR <- correlate(expression_AML_riboD_DR, method = "spearman", quiet = TRUE)
```

``` r
# ALL_polyA
sum_cor_ALL_polyA_DR <- cor_ALL_polyA_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_ALL_polyA_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_1667_S01 | 0.5901758 | 0.9434441 |    0.8793959 |  0.8694860 | 0.0457651 |
| THR24_2131_S01 | 0.5938943 | 0.9661436 |    0.8721808 |  0.8570567 | 0.0553884 |
| THR24_2119_S01 | 0.5751894 | 0.9534580 |    0.8693387 |  0.8573294 | 0.0547977 |
| THR24_1921_S01 | 0.5780204 | 0.9617942 |    0.8675398 |  0.8631504 | 0.0507880 |
| THR24_1845_S01 | 0.4949929 | 0.9758476 |    0.8672127 |  0.8564498 | 0.0616542 |

``` r
# ALL_riboD
sum_cor_ALL_riboD_DR <- cor_ALL_riboD_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_ALL_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_4203_S01 | 0.6106342 | 0.9701557 |    0.8821076 |  0.8729732 | 0.0496029 |
| THR24_3471_S01 | 0.6473365 | 0.9628727 |    0.8777846 |  0.8678079 | 0.0501633 |
| THR24_3688_S01 | 0.6198272 | 0.9701557 |    0.8759759 |  0.8686471 | 0.0529331 |
| THR24_4235_S01 | 0.6608003 | 0.9777067 |    0.8748054 |  0.8679768 | 0.0484338 |
| THR24_3828_S01 | 0.6060137 | 0.9500586 |    0.8744380 |  0.8669066 | 0.0488452 |

``` r
# AML_polyA
sum_cor_AML_polyA_DR <- cor_AML_polyA_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_AML_polyA_DR, n = 5)
```

| term            | min_value | max_value | median_value | mean_value |     stdev |
|:----------------|----------:|----------:|-------------:|-----------:|----------:|
| TCGA-AB-2889-03 | 0.5255989 | 0.9787678 |    0.9009624 |  0.8915278 | 0.0522497 |
| TCGA-AB-2844-03 | 0.5529009 | 0.9588590 |    0.9005687 |  0.8903770 | 0.0511095 |
| TCGA-AB-2846-03 | 0.5355347 | 0.9787678 |    0.9003138 |  0.8920143 | 0.0525085 |
| TCGA-AB-2881-03 | 0.5336866 | 0.9788672 |    0.9002967 |  0.8906179 | 0.0518842 |
| TCGA-AB-2996-03 | 0.5406430 | 0.9595678 |    0.8999107 |  0.8843799 | 0.0579246 |

``` r
# AML_riboD
sum_cor_AML_riboD_DR <- cor_AML_riboD_DR %>%
  pivot_longer(cols = -term, # convert to long format, reshapes the data from wide (50 columns) to long format (one row per sample pair)
               names_to = "other_sample",
               values_to = "correlation") %>%
  filter(term != other_sample) %>%  #removes the self-correlations
  group_by(term) %>%
  summarise(
    min_value = min(correlation, na.rm = TRUE),
    max_value = max(correlation, na.rm = TRUE),
    median_value = median(correlation, na.rm = TRUE),
    mean_value = mean(correlation, na.rm = TRUE),
    stdev = sd(correlation, na.rm = TRUE)
  ) %>%
  arrange(desc(median_value))
head(sum_cor_AML_riboD_DR, n = 5)
```

| term           | min_value | max_value | median_value | mean_value |     stdev |
|:---------------|----------:|----------:|-------------:|-----------:|----------:|
| THR24_4050_S01 | 0.5628071 | 0.9434897 |    0.8824564 |  0.8684954 | 0.0594492 |
| THR24_4356_S01 | 0.4495990 | 0.9728473 |    0.8809816 |  0.8666974 | 0.0730056 |
| THR24_4265_S01 | 0.6100284 | 0.9500586 |    0.8808372 |  0.8688932 | 0.0612926 |
| THR24_4263_S01 | 0.6061449 | 0.9500586 |    0.8775751 |  0.8687184 | 0.0573540 |
| THR24_4352_S01 | 0.5402336 | 0.9484302 |    0.8756370 |  0.8684385 | 0.0654422 |

``` r
filtered_ALL_polyA_list <- sum_cor_ALL_polyA_DR %>%
  slice_head(n = 20) %>%
  select(term) %>%
  mutate(disease_and_prep = "ALL_polyA")

filtered_ALL_riboD_list <- sum_cor_ALL_riboD_DR %>%
  slice_head(n = 20) %>% # taking the top most correlated samples
  select(term) %>%
  mutate(disease_and_prep = "ALL_riboD") # need to add back in the category

filtered_AML_polyA_list <- sum_cor_AML_polyA_DR %>%
  slice_head(n = 20) %>%
  select(term) %>%
  mutate(disease_and_prep = "AML_polyA")

filtered_AML_riboD_list <- sum_cor_AML_riboD_DR %>%
  slice_head(n = 20) %>%
  select(term) %>%
  mutate(disease_and_prep = "AML_riboD")
```

``` r
filtered_ALL_AML_list <- bind_rows(filtered_ALL_polyA_list, filtered_ALL_riboD_list, filtered_AML_polyA_list, filtered_AML_riboD_list)

write_tsv(filtered_ALL_AML_list, "../input_data/filtered_ALL_AML_list.tsv")
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
     date     2026-06-08
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     ! package      * version date (UTC) lib source
     P bit            4.6.0   2025-03-06 [?] RSPM
     P bit64          4.8.0   2026-04-21 [?] RSPM
     P cli            3.6.6   2026-04-09 [?] RSPM
     P corrr        * 0.4.5   2025-08-18 [?] RSPM
     P crayon         1.5.3   2024-06-20 [?] RSPM
     P curl           7.1.0   2026-04-22 [?] RSPM
     P digest         0.6.39  2025-11-19 [?] RSPM
     P dplyr        * 1.2.1   2026-04-03 [?] RSPM
     P evaluate       1.0.5   2025-08-27 [?] RSPM
     P farver         2.1.2   2024-05-13 [?] RSPM
     P fastmap        1.2.0   2024-05-15 [?] RSPM
     P forcats      * 1.0.1   2025-09-25 [?] RSPM
     P generics       0.1.4   2025-05-09 [?] RSPM
     P ggplot2      * 4.0.3   2026-04-22 [?] RSPM
     P glue           1.8.1   2026-04-17 [?] RSPM
     P gtable         0.3.6   2024-10-25 [?] RSPM
     P hms            1.1.4   2025-10-17 [?] RSPM
     P htmltools      0.5.9   2025-12-04 [?] RSPM
     P jsonlite       2.0.0   2025-03-27 [?] RSPM
     P knitr          1.51    2025-12-20 [?] RSPM
     P lifecycle      1.0.5   2026-01-08 [?] RSPM
     P lubridate    * 1.9.5   2026-02-04 [?] RSPM
     P magrittr       2.0.5   2026-04-04 [?] RSPM
     P pillar         1.11.1  2025-09-17 [?] RSPM
     P pkgconfig      2.0.3   2019-09-22 [?] RSPM
     P purrr        * 1.2.2   2026-04-10 [?] RSPM
     P R6             2.6.1   2025-02-15 [?] RSPM
     P RColorBrewer   1.1-3   2022-04-03 [?] RSPM
     P readr        * 2.2.0   2026-02-19 [?] RSPM
       renv           1.1.5   2025-07-24 [1] RSPM (R 4.5.2)
     P rlang          1.2.0   2026-04-06 [?] RSPM
     P rmarkdown      2.31    2026-03-26 [?] RSPM
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
     P xfun           0.57    2026-03-20 [?] RSPM
     P yaml           2.3.12  2025-12-10 [?] RSPM

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/correlation_analysis/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Users/maryke/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
