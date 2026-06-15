# Fig1G


## Fig 1G

Scatterplot showing median expression per gene.

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
library(ggrepel)
```

``` r
# sample ID files
# polyA
SS_polyA_list <- read_tsv("../../input_data/sample_selection/SS_polyA_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
aRMS_polyA_list <- read_tsv("../../input_data/sample_selection/aRMS_polyA_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_polyA_list <- read_tsv("../../input_data/sample_selection/WT_polyA_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NB_polyA_list <- read_tsv("../../input_data/sample_selection/NB_polyA_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_polyA_list <- read_tsv("../../input_data/sample_selection/ALL_polyA_list.tsv")
```

    Rows: 20 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AML_polyA_list <- read_tsv("../../input_data/sample_selection/AML_polyA_list.tsv")
```

    Rows: 20 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# riboD
SS_riboD_list <- read_tsv("../../input_data/sample_selection/SS_riboD_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
aRMS_riboD_list <- read_tsv("../../input_data/sample_selection/aRMS_riboD_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_riboD_list <- read_tsv("../../input_data/sample_selection/WT_riboD_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NB_riboD_list <- read_tsv("../../input_data/sample_selection/NB_riboD_list.tsv")
```

    Rows: 15 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_riboD_list <- read_tsv("../../input_data/sample_selection/ALL_riboD_list.tsv")
```

    Rows: 20 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AML_riboD_list <- read_tsv("../../input_data/sample_selection/AML_riboD_list.tsv")
```

    Rows: 20 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): term, disease_and_prep

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# log2(TPM+1) expression files
# polyA
SS_polyA_expr <- read_tsv("../../input_data/sample_selection/SS_polyA_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR39_1373_S01, TCGA-WK-A8XT-01, TH40_2281_S01, THR39_1375_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
aRMS_polyA_expr <- read_tsv("../../input_data/sample_selection/aRMS_polyA_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR29_0788_S01, THR29_0775_S01, THR29_0757_S01, THR29_0762_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_polyA_expr <- read_tsv("../../input_data/sample_selection/WT_polyA_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): TARGET-50-PAJNCZ-01, TARGET-50-PAEBXA-01, TARGET-50-PALERC-01, TAR...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NB_polyA_expr <- read_tsv("../../input_data/sample_selection/NB_polyA_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): TARGET-30-PASUML-01, TARGET-30-PASEGA-01, TARGET-30-PAPUAR-01, TAR...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_polyA_expr <- read_tsv("../../input_data/sample_selection/ALL_polyA_expr.tsv")
```

    Rows: 60498 Columns: 21
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (20): THR24_1667_S01, THR24_2131_S01, THR24_2119_S01, THR24_1921_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AML_polyA_expr <- read_tsv("../../input_data/sample_selection/AML_polyA_expr.tsv")
```

    Rows: 60498 Columns: 21
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (20): TCGA-AB-2889-03, TCGA-AB-2844-03, TCGA-AB-2846-03, TCGA-AB-2881-03...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# riboD
SS_riboD_expr <- read_tsv("../../input_data/sample_selection/SS_riboD_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR51_4556_S01, THR51_4558_S01, THR51_4554_S01, THR24_3992_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
aRMS_riboD_expr <- read_tsv("../../input_data/sample_selection/aRMS_riboD_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR24_3244_S01, THR24_3371_S01, THR24_3181_S01, THR24_3178_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
WT_riboD_expr <- read_tsv("../../input_data/sample_selection/WT_riboD_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR24_3218_S01, THR24_4194_S01, THR24_4284_S01, THR24_4369_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NB_riboD_expr <- read_tsv("../../input_data/sample_selection/NB_riboD_expr.tsv")
```

    Rows: 60498 Columns: 16
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (15): THR24_4310_S01, THR24_2779_S01, THR24_3516_S01, THR24_4114_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ALL_riboD_expr <- read_tsv("../../input_data/sample_selection/ALL_riboD_expr.tsv")
```

    Rows: 60498 Columns: 21
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (20): THR24_4203_S01, THR24_3471_S01, THR24_3688_S01, THR24_4235_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
AML_riboD_expr <- read_tsv("../../input_data/sample_selection/AML_riboD_expr.tsv")
```

    Rows: 60498 Columns: 21
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (1): Gene
    dbl (20): THR24_4050_S01, THR24_4356_S01, THR24_4265_S01, THR24_4263_S01, TH...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# to convert EnsemblIDs to HugoIDs
gene_names <- read.table("../../input_data/EnsGeneID_Hugo_Observed_Conversions.txt",
header = TRUE, sep = "\t", stringsAsFactors = FALSE )
```

Compute gene medians

``` r
compute_gene_medians <- function(polyA_df, riboD_df, disease_name) {
  polyA_Medians <- polyA_df %>%
    rowwise() %>%
    mutate(polyA_medians = median(c_across(-Gene), na.rm = TRUE)) %>%
    ungroup() %>%
    select(Gene, polyA_medians)

  riboD_Medians <- riboD_df %>%
    rowwise() %>%
    mutate(riboD_medians = median(c_across(-Gene), na.rm = TRUE)) %>%
    ungroup() %>%
    select(Gene, riboD_medians)

  merged <- polyA_Medians %>%
    inner_join(riboD_Medians, by = "Gene") %>%
    mutate(Disease = disease_name)

  return(merged)
}
```

Applying function to each disease/library prep type

``` r
aRMS_medians <- compute_gene_medians(aRMS_polyA_expr, aRMS_riboD_expr, "aRMS")
SS_medians   <- compute_gene_medians(SS_polyA_expr,   SS_riboD_expr,   "SS")
WT_medians   <- compute_gene_medians(WT_polyA_expr,   WT_riboD_expr,   "WT")
NB_medians   <- compute_gene_medians(NB_polyA_expr,   NB_riboD_expr,   "NB")
ALL_medians   <- compute_gene_medians(ALL_polyA_expr,   ALL_riboD_expr,   "ALL")
AML_medians   <- compute_gene_medians(AML_polyA_expr,   AML_riboD_expr,   "AML")

combined_medians <- bind_rows(aRMS_medians, SS_medians, WT_medians, NB_medians, ALL_medians, AML_medians)
```

``` r
plot_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      plot.title = element_text(face = "bold", size = base_size * 1, hjust = 0.5),
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

convert EnsemblIDs to Hugo IDs for readability

``` r
convert_to_hugo <- function(df) {
  df %>%
    left_join(gene_names, by = c("Gene" = "EnsGeneID")) %>%
    mutate(Gene = ifelse(!is.na(HugoID), HugoID, Gene)) %>%
    select(-HugoID)
}

combined_medians_hugo <- combined_medians %>%
  convert_to_hugo()
```

Compute residual for each disease

``` r
combined_medians_hugo_res <- combined_medians_hugo %>%
  mutate(
    residual = riboD_medians - polyA_medians
  )
```

Define representative polyadenylated and non-polyadenylated genes

``` r
# representative polyadenylated genes were found in Fig1E_F, calculated as the median expression level in PolyA samples divided by the median expression level in RiboD samples


poly_genes <- c(
  SS  = "MTMR2",
  aRMS = "TM2D1",
  WT   = "SYPL1",
  NB   = "CCDC64",
  ALL     = "XPO5",
  AML     = "EML3"
    
)
nonpoly_gene <- "HIST1H1B"
```

Color labeled genes

``` r
combined_medians_hugo_res <- combined_medians_hugo_res %>%
  mutate(
    label_color = case_when(
      Gene %in% poly_genes[Disease] ~ "#BB5566",                     # pinkish red
      Gene %in% nonpoly_gene ~ "#BB5566",                         # pinkish red
      Gene == poly_genes & Disease == "ALL" ~ "#BB5566",
      TRUE ~ "grey70"
    )
  )
```

Create labeled dateframe with representative genes

``` r
combined_medians_hugo_res <- combined_medians_hugo_res %>%
  mutate(
    label = case_when(
      Gene == poly_genes[Disease] ~ Gene,
      Gene == nonpoly_gene ~ Gene,            # HIST1H1B + median ratio gene
      TRUE ~ NA_character_
    )
  )
```

Plot average expression across library prep methods for each disease

``` r
diseases <- unique(combined_medians_hugo_res$Disease)

scatterplots_by_disease <- map(diseases, function(d) {

  df <- combined_medians_hugo_res %>% filter(Disease == d)

 ggplot(df, aes(x = polyA_medians, y = riboD_medians)) +
    
   # background points
    geom_point(data = df %>% filter(is.na(label)),
               color = "grey70", alpha = 0.5, size = 1.2) +

    # highlighted labeled points
    geom_point(data = df %>% filter(!is.na(label)),
               aes(color = label_color),
               size = 2.5, alpha = 0.9) +

    
       # Labels for representative genes and HIST1H1B
    geom_text_repel(
      aes(label = label),
      size = 4.5,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.4,
      point.padding = 0.2
    ) +

    scale_color_identity() +   # use colors exactly as defined

    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "gray50") +

    labs(
      title = paste(d),
      x = "Average PolyA Expression (log2(TPM+1))",
      y = "Average RiboD Expression (log2(TPM+1))",
    ) +
    theme_minimal(base_size = 14) +
    plot_theme()
})
names(scatterplots_by_disease) <- diseases

scatterplots_by_disease
```

    $aRMS

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-1.png)


    $SS

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-2.png)


    $WT

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-3.png)


    $NB

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-4.png)


    $ALL

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-5.png)


    $AML

    Warning: Removed 60496 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1G_files/figure-commonmark/Fig1G-6.png)
