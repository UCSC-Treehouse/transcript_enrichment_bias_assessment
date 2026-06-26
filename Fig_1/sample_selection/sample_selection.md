# sample_selection


## Selecting samples based on correlation analysis

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
library(janitor)
```


    Attaching package: 'janitor'

    The following objects are masked from 'package:stats':

        chisq.test, fisher.test

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
SS_aRMS_counts <- read_tsv("../../input_data/expression_SS_aRMS.tsv.gz")
```

    Rows: 60498 Columns: 170
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr   (1): Gene
    dbl (169): TH03_0112_S01, TH03_0112_S02, TH27_1157_S01, TH27_1157_S02, TH27_...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_NB_counts <- read_tsv("../../input_data/expression_WT_NB.tsv.gz")
```

    Rows: 60498 Columns: 502
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr   (1): Gene
    dbl (501): TH03_0012_S01, TH03_0012_S02, TH27_1356_S01, TH27_1403_S01, TH27_...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_AML_counts <- read_tsv("../../input_data/expression_ALL_AML.tsv.gz")
```

    Rows: 60498 Columns: 1927
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr    (1): Gene
    dbl (1926): TH03_0011_S01, TH03_0118_S01, TH34_1238_S01, TH34_1239_S01, TH27...

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
# SS polyA
SS_polyA_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("SS_polyA")) %>% 
  select(term, disease_and_prep)

SS_polyA_log2tpm1 <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(SS_polyA_list$term))

SS_polyA_counts <- SS_aRMS_counts %>%
  select(Gene, all_of(SS_polyA_list$term))

# SS riboD
SS_riboD_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("SS_riboD")) %>% 
  select(term, disease_and_prep)

SS_riboD_log2tpm1 <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(SS_riboD_list$term))

SS_riboD_counts <- SS_aRMS_counts %>%
  select(Gene, all_of(SS_riboD_list$term))

write_tsv(SS_polyA_list, "../../input_data/sample_selection/SS_polyA_list.tsv")
write_tsv(SS_riboD_list, "../../input_data/sample_selection/SS_riboD_list.tsv")
write_tsv(SS_polyA_log2tpm1, "../../input_data/sample_selection/SS_polyA_log2tpm1.tsv")
write_tsv(SS_riboD_log2tpm1, "../../input_data/sample_selection/SS_riboD_log2tpm1.tsv")
write_tsv(SS_polyA_counts, "../../input_data/sample_selection/SS_polyA_counts.tsv")
write_tsv(SS_riboD_counts, "../../input_data/sample_selection/SS_riboD_counts.tsv")

# aRMS polyA
aRMS_polyA_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("aRMS_polyA")) %>% 
  select(term, disease_and_prep)

aRMS_polyA_log2tpm1 <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(aRMS_polyA_list$term))

aRMS_polyA_counts <- SS_aRMS_counts %>%
  select(Gene, all_of(aRMS_polyA_list$term))

# aRMS riboD
aRMS_riboD_list <- filtered_SS_aRMS_list %>%
  filter(disease_and_prep %in% c("aRMS_riboD")) %>% 
  select(term, disease_and_prep)

aRMS_riboD_log2tpm1 <- SS_aRMS_log2tpm1 %>%
  select(Gene, all_of(aRMS_riboD_list$term))

aRMS_riboD_counts <- SS_aRMS_counts %>%
  select(Gene, all_of(aRMS_riboD_list$term))

write_tsv(aRMS_polyA_list, "../../input_data/sample_selection/aRMS_polyA_list.tsv")
write_tsv(aRMS_riboD_list, "../../input_data/sample_selection/aRMS_riboD_list.tsv")
write_tsv(aRMS_polyA_log2tpm1, "../../input_data/sample_selection/aRMS_polyA_log2tpm1.tsv")
write_tsv(aRMS_riboD_log2tpm1, "../../input_data/sample_selection/aRMS_riboD_log2tpm1.tsv")
write_tsv(aRMS_polyA_counts, "../../input_data/sample_selection/aRMS_polyA_counts.tsv")
write_tsv(aRMS_riboD_counts, "../../input_data/sample_selection/aRMS_riboD_counts.tsv")

# WT polyA
WT_polyA_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("WT_polyA")) %>% 
  select(term, disease_and_prep)

WT_polyA_log2tpm1 <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(WT_polyA_list$term))

WT_polyA_counts <- WT_NB_counts %>%
  select(Gene, all_of(WT_polyA_list$term))

# WT riboD 
WT_riboD_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("WT_riboD")) %>% 
  select(term, disease_and_prep)

WT_riboD_log2tpm1 <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(WT_riboD_list$term))

WT_riboD_counts <- WT_NB_counts %>%
  select(Gene, all_of(WT_riboD_list$term))

write_tsv(WT_polyA_list, "../../input_data/sample_selection/WT_polyA_list.tsv")
write_tsv(WT_riboD_list, "../../input_data/sample_selection/WT_riboD_list.tsv")
write_tsv(WT_polyA_log2tpm1, "../../input_data/sample_selection/WT_polyA_log2tpm1.tsv")
write_tsv(WT_riboD_log2tpm1, "../../input_data/sample_selection/WT_riboD_log2tpm1.tsv")
write_tsv(WT_polyA_counts, "../../input_data/sample_selection/WT_polyA_counts.tsv")
write_tsv(WT_riboD_counts, "../../input_data/sample_selection/WT_riboD_counts.tsv")


# NB polyA
NB_polyA_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("NB_polyA")) %>% 
  select(term, disease_and_prep)

NB_polyA_log2tpm1 <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(NB_polyA_list$term))

NB_polyA_counts <- WT_NB_counts %>%
  select(Gene, all_of(NB_polyA_list$term))

# NB riboD
NB_riboD_list <- filtered_WT_NB_list %>%
  filter(disease_and_prep %in% c("NB_riboD")) %>% 
  select(term, disease_and_prep)

NB_riboD_log2tpm1 <- WT_NB_log2tpm1 %>%
  select(Gene, all_of(NB_riboD_list$term))

NB_riboD_counts <- WT_NB_counts %>%
  select(Gene, all_of(NB_riboD_list$term))

write_tsv(NB_polyA_list, "../../input_data/sample_selection/NB_polyA_list.tsv")
write_tsv(NB_riboD_list, "../../input_data/sample_selection/NB_riboD_list.tsv")
write_tsv(NB_polyA_log2tpm1, "../../input_data/sample_selection/NB_polyA_log2tpm1.tsv")
write_tsv(NB_riboD_log2tpm1, "../../input_data/sample_selection/NB_riboD_log2tpm1.tsv")
write_tsv(NB_polyA_counts, "../../input_data/sample_selection/NB_polyA_counts.tsv")
write_tsv(NB_riboD_counts, "../../input_data/sample_selection/NB_riboD_counts.tsv")

# ALL polyA
ALL_polyA_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("ALL_polyA")) %>% 
  select(term, disease_and_prep)

ALL_polyA_log2tpm1 <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(ALL_polyA_list$term))

ALL_polyA_counts <- ALL_AML_counts %>%
  select(Gene, all_of(ALL_polyA_list$term))

# ALL riboD
ALL_riboD_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("ALL_riboD")) %>% 
  select(term, disease_and_prep)

ALL_riboD_log2tpm1 <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(ALL_riboD_list$term))

ALL_riboD_counts <- ALL_AML_counts %>%
  select(Gene, all_of(ALL_riboD_list$term))

write_tsv(ALL_polyA_list, "../../input_data/sample_selection/ALL_polyA_list.tsv")
write_tsv(ALL_riboD_list, "../../input_data/sample_selection/ALL_riboD_list.tsv")
write_tsv(ALL_polyA_log2tpm1, "../../input_data/sample_selection/ALL_polyA_log2tpm1.tsv")
write_tsv(ALL_riboD_log2tpm1, "../../input_data/sample_selection/ALL_riboD_log2tpm1.tsv")
write_tsv(ALL_polyA_counts, "../../input_data/sample_selection/ALL_polyA_counts.tsv")
write_tsv(ALL_riboD_counts, "../../input_data/sample_selection/ALL_riboD_counts.tsv")

# AML polyA
AML_polyA_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("AML_polyA")) %>% 
  select(term, disease_and_prep)

AML_polyA_log2tpm1 <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(AML_polyA_list$term))

AML_polyA_counts <- ALL_AML_counts %>%
  select(Gene, all_of(AML_polyA_list$term))

# AML riboD
AML_riboD_list <- filtered_ALL_AML_list %>%
  filter(disease_and_prep %in% c("AML_riboD")) %>% 
  select(term, disease_and_prep)

AML_riboD_log2tpm1 <- ALL_AML_log2tpm1 %>%
  select(Gene, all_of(AML_riboD_list$term))

AML_riboD_counts <- ALL_AML_counts %>%
  select(Gene, all_of(AML_riboD_list$term))

write_tsv(AML_polyA_list, "../../input_data/sample_selection/AML_polyA_list.tsv")
write_tsv(AML_riboD_list, "../../input_data/sample_selection/AML_riboD_list.tsv")
write_tsv(AML_polyA_log2tpm1, "../../input_data/sample_selection/AML_polyA_log2tpm1.tsv")
write_tsv(AML_riboD_log2tpm1, "../../input_data/sample_selection/AML_riboD_log2tpm1.tsv")
write_tsv(AML_polyA_counts, "../../input_data/sample_selection/AML_polyA_counts.tsv")
write_tsv(AML_riboD_counts, "../../input_data/sample_selection/AML_riboD_counts.tsv")
```

### Samples for UMAP

``` r
subset_samples <- read_tsv("../../../lib-prep-visualization/matched_subsamples/polyA_riboD_disease_matched/subset_samples.tsv") %>% # this contains the samples that were randomly subsetted for UMAP
  add_row(THR19_0450_S01 = "THR19_0450_S01", .before = 1) %>% # subset_samples had THR19_0450_S01 as column name, not as data point
  rename(th_dataset_id = THR19_0450_S01)
```

    Rows: 3245 Columns: 1
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (1): THR19_0450_S01

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
polyA_v25 <- read_tsv("../../input_data/clinical_Treehouse-Tumor-Compendium-25.01-PolyA_20250131v1.tsv") %>% # will need this to determine which THIDs we can remove safely 
  mutate(compendia = "polyA") %>%
  select(th_dataset_id, disease, compendia)
```

    Rows: 13359 Columns: 12
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (12): th_dataset_id, age_at_dx, disease, icd_disease, organism, pedaya, ...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
riboD_v25 <- read_tsv("../../input_data/clinical_Treehouse-Tumor-Compendium-25.01-RiboD_20250306v1.tsv") %>%
  mutate(compendia = "riboD") %>%
  select(th_dataset_id, disease, compendia)
```

    Rows: 2079 Columns: 12
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (12): th_dataset_id, age_at_dx, disease, icd_disease, organism, pedaya, ...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# During our analysis, we discovered that some ALL samples in the PolyA compendium were actually derived from a ribo-depleted RNA-seq library and were mislabeled as PolyA. We remove those mislabeled ALL samples from the UMAP.
mislabeled_ALLpolyA <- read_tsv("../../input_data/TARGET_ALL_SSTR_SraRun.tsv") %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) %>% # this contains the THIDs that we need to remove from the PolyA samples, since they came from a riboD library. Many THIDs are duplicates.
  select(SAMPLE_ID, disease_and_prep) %>%
  rename(th_dataset_id = SAMPLE_ID)
```

    Rows: 118 Columns: 132
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (45): SAMPLE_ID, disease_and_prep, SUBJECT_ID, SEX, consent_abbreviatio...
    dbl   (9): CONSENT, dbgap_subject_id, dbgap_sample_id, Bytes, Consent_Code, ...
    lgl  (76): subject_is_affected, study_disease, alignment_software (exp), Com...
    dttm  (2): ReleaseDate, create_date

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# combining clinical files and samples we need to label
polyA_riboD <- rbind(polyA_v25, riboD_v25)

all_interested_samples <- rbind(filtered_SS_aRMS_list, filtered_WT_NB_list, filtered_ALL_AML_list) %>%
  rename(th_dataset_id = term)

write_tsv(all_interested_samples, "../../input_data/sample_selection/forUMAP/all_interested_samples.tsv")
```

``` r
subset_samples_disease <- left_join(subset_samples, polyA_riboD, by = "th_dataset_id")
# joining so I know which THIDs in the subset_samples are from polyA, which from riboD, which are from which disease
```

``` r
# how many mislabeled ALL_polyA samples were in the subset_samples?
subset_samples_disease_mislabeledALL <- subset_samples_disease %>%
  filter(subset_samples_disease$th_dataset_id %in% mislabeled_ALLpolyA$th_dataset_id)

# remove the mislabeled ALL_polyA samples from subset_samples_disease
subset_samples_disease_wo_mislabeledALL <- subset_samples_disease %>%
  filter(!subset_samples_disease$th_dataset_id %in% mislabeled_ALLpolyA$th_dataset_id)

print(paste("the number of samples in subset_samples is ",nrow(subset_samples_disease)))
```

    [1] "the number of samples in subset_samples is  3246"

``` r
print(paste("the number of mislabeled ALL_polyA samples in subset_samples is ", nrow(subset_samples_disease_mislabeledALL)))
```

    [1] "the number of mislabeled ALL_polyA samples in subset_samples is  22"

``` r
print(paste("the new number of samples in subset_samples after removing mislabeled is ", nrow(subset_samples_disease_wo_mislabeledALL)))
```

    [1] "the new number of samples in subset_samples after removing mislabeled is  3224"

``` r
# how many samples are in subset_samples (after removing mislabeled ALL) for each disease?
counts_woMislabeled <- subset_samples_disease_wo_mislabeledALL %>%
  group_by(disease, compendia) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
```

    `summarise()` has grouped output by 'disease'. You can override using the
    `.groups` argument.

``` r
head(counts_woMislabeled, 8)
```

| disease                      | compendia | n_disease |
|:-----------------------------|:----------|----------:|
| glioma                       | polyA     |       437 |
| glioma                       | riboD     |       437 |
| acute lymphoblastic leukemia | riboD     |       328 |
| acute lymphoblastic leukemia | polyA     |       306 |
| medulloblastoma              | polyA     |       125 |
| medulloblastoma              | riboD     |       125 |
| ependymoma                   | polyA     |       102 |
| ependymoma                   | riboD     |       102 |

``` r
# now we need to remove the same number of ALL_riboD samples (22) from subset_samples to balance it out
# but first we need to pull out the subset_samples samples that we are not interested in labelling

# anti_join returns all rows from the first data frame (x) that do not have a match in the second data frame

# the first dataframe will be all the ALL riboD samples from subset_samples_disease_wo_mislabeledALL
# the second dataframe will be ALL_riboD that we want to label

ALL_riboD_subset <- subset_samples_disease_wo_mislabeledALL %>%
  filter(disease == "acute lymphoblastic leukemia") %>%
  filter(compendia == "riboD")

ALL_riboD_tolabel <- filtered_ALL_AML_list %>%
  filter(disease_and_prep == "ALL_riboD") %>%
  rename(th_dataset_id = term)

ALL_riboD_can_remove <- anti_join(ALL_riboD_subset, ALL_riboD_tolabel, by = "th_dataset_id")

# finding the 22 ALL riboD samples we can remove
set.seed(123)
ALL_riboD_to_remove <- ALL_riboD_can_remove %>%
  sample_n(22)

# removing 22 ALL riboD samples from subset_samples
subset_samples_wo_polyAriboD_ALL <- subset_samples_disease_wo_mislabeledALL %>%
  filter(!subset_samples_disease_wo_mislabeledALL$th_dataset_id %in% ALL_riboD_to_remove$th_dataset_id)

print(paste("the number of samples in subset_samples after removing 22 mislabeled ALL_polyA is ",nrow(subset_samples_disease_wo_mislabeledALL)))
```

    [1] "the number of samples in subset_samples after removing 22 mislabeled ALL_polyA is  3224"

``` r
print(paste("the number of ALL riboD samples we want to label is ",nrow(ALL_riboD_tolabel)))
```

    [1] "the number of ALL riboD samples we want to label is  20"

``` r
print(paste("the number of ALL riboD samples that was present in subset is ",nrow(ALL_riboD_subset)))
```

    [1] "the number of ALL riboD samples that was present in subset is  328"

``` r
print(paste("the number of ALL riboD samples we are not interested in labelling is ",nrow(ALL_riboD_can_remove)))
```

    [1] "the number of ALL riboD samples we are not interested in labelling is  308"

``` r
print(paste("the number of subset_samples after removing the 22 mislabeled ALL_polyA samples and corresponding 22 ALL_riboD samples is ",nrow(subset_samples_wo_polyAriboD_ALL)))
```

    [1] "the number of subset_samples after removing the 22 mislabeled ALL_polyA samples and corresponding 22 ALL_riboD samples is  3202"

``` r
# double check that the sample counts are the same for each disease after removing 22 ALL samples from polyA and riboD

counts_subset_samples_wo_polyAriboD_ALL <- subset_samples_wo_polyAriboD_ALL %>%
  group_by(disease, compendia) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
```

    `summarise()` has grouped output by 'disease'. You can override using the
    `.groups` argument.

``` r
head(counts_subset_samples_wo_polyAriboD_ALL, 8)
```

| disease                      | compendia | n_disease |
|:-----------------------------|:----------|----------:|
| glioma                       | polyA     |       437 |
| glioma                       | riboD     |       437 |
| acute lymphoblastic leukemia | polyA     |       306 |
| acute lymphoblastic leukemia | riboD     |       306 |
| medulloblastoma              | polyA     |       125 |
| medulloblastoma              | riboD     |       125 |
| ependymoma                   | polyA     |       102 |
| ependymoma                   | riboD     |       102 |

``` r
# anti_join returns all rows from the first data frame (x) that do not have a match in the second data frame
samples_to_add <- anti_join(all_interested_samples, subset_samples_wo_polyAriboD_ALL, by = "th_dataset_id") # these are the samples that we're interested in labeling but are not represented in subset_samples (after removing ALL)


# how many of each disease do we need to add?
samples_to_add_counts <- samples_to_add %>%
  group_by(disease_and_prep) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
samples_to_add_counts
```

| disease_and_prep | n_disease |
|:-----------------|----------:|
| AML_polyA        |        17 |
| ALL_polyA        |        15 |
| aRMS_polyA       |        13 |
| WT_polyA         |        11 |
| NB_polyA         |         8 |
| SS_polyA         |         6 |

``` r
samples_can_remove <- anti_join(subset_samples_wo_polyAriboD_ALL, all_interested_samples, by = "th_dataset_id") %>%
  filter(compendia == "polyA") # these are samples that we can safely remove from subset_samples after adding the samples we need (in order to balance the sample numbers)

# anti_join returns all rows from the first data frame (x) that do not have a match in the second data frame

samples_can_remove_counts <- samples_can_remove %>%
  group_by(disease, compendia) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
```

    `summarise()` has grouped output by 'disease'. You can override using the
    `.groups` argument.

``` r
head(samples_can_remove_counts, 8)
```

| disease                      | compendia | n_disease |
|:-----------------------------|:----------|----------:|
| glioma                       | polyA     |       437 |
| acute lymphoblastic leukemia | polyA     |       301 |
| medulloblastoma              | polyA     |       125 |
| ependymoma                   | polyA     |       102 |
| acute myeloid leukemia       | polyA     |        96 |
| neuroblastoma                | polyA     |        73 |
| wilms tumor                  | polyA     |        52 |
| osteosarcoma                 | polyA     |        49 |

``` r
# make a list of equal number of samples to remove from subset_samples
set.seed(1100)
ALL_toremove <- samples_can_remove %>%
  filter(disease == 'acute lymphoblastic leukemia') %>%
  sample_n(15)
  
AML_toremove <- samples_can_remove %>% 
  filter(disease == 'acute myeloid leukemia') %>%
  sample_n(17)

NB_toremove <- samples_can_remove %>% 
  filter(disease == 'neuroblastoma') %>% 
  sample_n(8)

SS_toremove <- samples_can_remove %>%
  filter(disease == 'synovial sarcoma') %>% 
  sample_n(6)

WT_toremove <- samples_can_remove %>%
  filter(disease == 'wilms tumor') %>%
  sample_n(11)

aRMS_toremove <- samples_can_remove %>%
  filter(disease == 'alveolar rhabdomyosarcoma') %>%
  sample_n(13)

samples_to_remove <- rbind(ALL_toremove, AML_toremove, NB_toremove, SS_toremove, WT_toremove, aRMS_toremove)
```

By now, we have removed the mislabeled ALL polyA samples and equal
number ALL_riboD samples from subset_samples, we have made a list of
samples we need to label that are missing from subset, and we have made
a list of samples that we can remove from subset after adding the
samples we want to label.

Next, we need to add the samples that we want to label that are not in
subset_samples, and then remove the same number of samples from
subset_samples to balance the numbers (or vice versa, we can remove
first and then add).

``` r
# removing the samples from each disease
subset_samples_wo_polyAriboD_ALL_removed <- subset_samples_wo_polyAriboD_ALL %>%
  filter(!subset_samples_wo_polyAriboD_ALL$th_dataset_id %in% samples_to_remove$th_dataset_id)

# double check the numbers - polyA should have n less samples than riboD for the specified diseases
subset_samples_wo_polyAriboD_ALL_removed_counts <- subset_samples_wo_polyAriboD_ALL_removed %>%
  group_by(disease, compendia) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
```

    `summarise()` has grouped output by 'disease'. You can override using the
    `.groups` argument.

``` r
head(subset_samples_wo_polyAriboD_ALL_removed_counts, 8)
```

| disease                      | compendia | n_disease |
|:-----------------------------|:----------|----------:|
| glioma                       | polyA     |       437 |
| glioma                       | riboD     |       437 |
| acute lymphoblastic leukemia | riboD     |       306 |
| acute lymphoblastic leukemia | polyA     |       291 |
| medulloblastoma              | polyA     |       125 |
| medulloblastoma              | riboD     |       125 |
| ependymoma                   | polyA     |       102 |
| ependymoma                   | riboD     |       102 |

``` r
# first, we need to reformat the samples_to_add because it currently has disease_and_prep and we need it in disease and compendia

samples_to_add_format <- samples_to_add %>%
  separate(col = disease_and_prep, into = c("disease", "compendia"), sep = "_") %>%
  mutate(disease = case_when(
    disease == "SS" ~ "synovial sarcoma",
    disease == "aRMS" ~ "alveolar rhabdomyosarcoma",
    disease == "WT" ~ "wilms tumor",
    disease == "NB" ~ "neuroblastoma",
    disease == "ALL" ~ "acute lymphoblastic leukemia",
    disease == "AML" ~ "acute myeloid leukemia"
  ))
```

``` r
# adding back the samples we want to label
subset_samples_new <- rbind(subset_samples_wo_polyAriboD_ALL_removed, samples_to_add_format)
```

``` r
# double check the numbers - polyA and riboD should have the same numbers now
subset_samples_new_counts <- subset_samples_new %>%
  group_by(disease, compendia) %>%
  summarise(n_disease = n()) %>% 
  arrange(desc(n_disease))
```

    `summarise()` has grouped output by 'disease'. You can override using the
    `.groups` argument.

``` r
head(subset_samples_new_counts, 8)
```

| disease                      | compendia | n_disease |
|:-----------------------------|:----------|----------:|
| glioma                       | polyA     |       437 |
| glioma                       | riboD     |       437 |
| acute lymphoblastic leukemia | polyA     |       306 |
| acute lymphoblastic leukemia | riboD     |       306 |
| medulloblastoma              | polyA     |       125 |
| medulloblastoma              | riboD     |       125 |
| ependymoma                   | polyA     |       102 |
| ependymoma                   | riboD     |       102 |

``` r
# saving the file for UMAP
write_tsv(subset_samples_new, "../../input_data/sample_selection/forUMAP/final_subset_samples_forUMAP.tsv")

subset_samples_new_IDonly <- subset_samples_new %>%
  select(th_dataset_id) %>%
  row_to_names(row_number = 1)

write_tsv(subset_samples_new_IDonly, "../../input_data/sample_selection/forUMAP/final_subset_samples_forUMAP_IDonly.tsv")
```

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
     date     2026-06-25
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     package      * version date (UTC) lib source
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
     bit64          4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
     cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
     digest         0.6.37  2024-08-19 [1] CRAN (R 4.5.0)
     dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
     forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.5.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
     ggplot2      * 4.0.0   2025-09-11 [1] CRAN (R 4.5.0)
     glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.5.0)
     htmltools      0.5.8.1 2024-04-04 [1] CRAN (R 4.5.0)
     janitor      * 2.2.1   2024-12-22 [1] CRAN (R 4.5.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
     knitr          1.50    2025-03-16 [1] CRAN (R 4.5.0)
     lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
     lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
     magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
     purrr        * 1.1.0   2025-07-10 [1] CRAN (R 4.5.0)
     R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
     RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
     readr        * 2.1.5   2024-01-10 [1] CRAN (R 4.5.0)
     rlang          1.2.0   2026-04-06 [1] CRAN (R 4.5.2)
     rmarkdown      2.30    2025-09-28 [1] CRAN (R 4.5.0)
     rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
     S7             0.2.0   2024-11-07 [1] CRAN (R 4.5.0)
     scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
     sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
     snakecase      0.11.1  2023-08-27 [1] CRAN (R 4.5.0)
     stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
     stringr      * 1.5.2   2025-09-08 [1] CRAN (R 4.5.0)
     tibble       * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
     tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.0)
     tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
     tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
     timechange     0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
     tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
     vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
     vroom          1.6.6   2025-09-19 [1] CRAN (R 4.5.0)
     withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
     xfun           0.55    2025-12-16 [1] CRAN (R 4.5.2)
     yaml           2.3.10  2024-07-26 [1] CRAN (R 4.5.0)

     [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
