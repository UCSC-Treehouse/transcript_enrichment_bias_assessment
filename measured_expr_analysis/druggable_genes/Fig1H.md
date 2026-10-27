# Fig1H


- [Fig 1H](#fig-1h)

## Fig 1H

Median expression of Treehouse Druggable Genes.

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
library(patchwork)
library(cowplot)
```


    Attaching package: 'cowplot'

    The following object is masked from 'package:patchwork':

        align_plots

    The following object is masked from 'package:lubridate':

        stamp

``` r
# sample ID files
# polyA
SS_polyA_list <- read_tsv("../../input_data/sample_selection/SS_polyA_list.tsv", show_col_types = FALSE)
aRMS_polyA_list <- read_tsv("../../input_data/sample_selection/aRMS_polyA_list.tsv", show_col_types = FALSE)
WT_polyA_list <- read_tsv("../../input_data/sample_selection/WT_polyA_list.tsv", show_col_types = FALSE)
NB_polyA_list <- read_tsv("../../input_data/sample_selection/NB_polyA_list.tsv", show_col_types = FALSE)
ALL_polyA_list <- read_tsv("../../input_data/sample_selection/ALL_polyA_list.tsv", show_col_types = FALSE)
AML_polyA_list <- read_tsv("../../input_data/sample_selection/AML_polyA_list.tsv", show_col_types = FALSE)

# riboD
SS_riboD_list <- read_tsv("../../input_data/sample_selection/SS_riboD_list.tsv", show_col_types = FALSE)
aRMS_riboD_list <- read_tsv("../../input_data/sample_selection/aRMS_riboD_list.tsv", show_col_types = FALSE)
WT_riboD_list <- read_tsv("../../input_data/sample_selection/WT_riboD_list.tsv", show_col_types = FALSE)
NB_riboD_list <- read_tsv("../../input_data/sample_selection/NB_riboD_list.tsv", show_col_types = FALSE)
ALL_riboD_list <- read_tsv("../../input_data/sample_selection/ALL_riboD_list.tsv", show_col_types = FALSE)
AML_riboD_list <- read_tsv("../../input_data/sample_selection/AML_riboD_list.tsv", show_col_types = FALSE)

# log2(TPM+1) expression files
SS_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/SS_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

SS_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/SS_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

aRMS_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/aRMS_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

aRMS_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/aRMS_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

WT_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/WT_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)
WT_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/WT_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

NB_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/NB_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)
NB_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/NB_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

ALL_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/ALL_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)
ALL_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/ALL_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)

AML_polyA_log2tpm1 <- read_tsv("../../input_data/sample_selection/AML_polyA_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)
AML_riboD_log2tpm1 <- read_tsv("../../input_data/sample_selection/AML_riboD_log2tpm1.tsv", show_col_types = FALSE) %>%
  rename(EnsGeneID = Gene)


# to convert EnsemblIDs to HugoIDs
gene_names <- read.table("../../input_data/EnsGeneID_Hugo_Observed_Conversions.txt",
header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# treehouse druggable gene list
drug_genes <- read_tsv("../../input_data/treehouseDruggableGenes_2019-06-12.txt", show_col_types = FALSE) %>%
  select(gene)
```

``` r
# getting the EnsemblIDs for the Treehouse druggable gene list
drug_gene_names <- drug_genes %>%
  rename(HugoID = gene) %>%
  left_join(gene_names, by = "HugoID")
```

Function for filtering for druggable genes

``` r
filter_druggable <- function(df, drug_gene_names) {
  df %>% filter(EnsGeneID %in% drug_gene_names$EnsGeneID)
    # left_join(drug_gene_names, by = "EnsGeneID") %>%
    # relocate(HugoID)
}
```

``` r
# filtering for TH druggable genes
aRMS_polyA_drug <- filter_druggable(aRMS_polyA_log2tpm1, drug_gene_names)
aRMS_riboD_drug <- filter_druggable(aRMS_riboD_log2tpm1, drug_gene_names)

SS_polyA_drug <- filter_druggable(SS_polyA_log2tpm1, drug_gene_names)
SS_riboD_drug <- filter_druggable(SS_riboD_log2tpm1, drug_gene_names)

WT_polyA_drug <- filter_druggable(WT_polyA_log2tpm1, drug_gene_names)
WT_riboD_drug <- filter_druggable(WT_riboD_log2tpm1, drug_gene_names)

NB_polyA_drug <- filter_druggable(NB_polyA_log2tpm1, drug_gene_names)
NB_riboD_drug <- filter_druggable(NB_riboD_log2tpm1, drug_gene_names)

ALL_polyA_drug <- filter_druggable(ALL_polyA_log2tpm1, drug_gene_names)
ALL_riboD_drug <- filter_druggable(ALL_riboD_log2tpm1, drug_gene_names)

AML_polyA_drug <- filter_druggable(AML_polyA_log2tpm1, drug_gene_names)
AML_riboD_drug <- filter_druggable(AML_riboD_log2tpm1, drug_gene_names)
```

Function for calculating median expression across samples

``` r
compute_gene_medians_polyA <- function(polyA_df, disease_name, Compendia ) {
  polyA_medians <- polyA_df %>%
    rowwise() %>%
    mutate(Expression = median(c_across(-EnsGeneID), na.rm = TRUE)) %>%
    ungroup() %>%
    select(EnsGeneID, Expression) %>%
    left_join(drug_gene_names, by = "EnsGeneID") %>%
    relocate(HugoID)
}

compute_gene_medians_riboD <- function(riboD_df, disease_name, Compendia) {
   riboD_medians <- riboD_df %>%
    rowwise() %>%
    mutate(Expression = median(c_across(-EnsGeneID), na.rm = TRUE)) %>%
    ungroup() %>%
    select(EnsGeneID, Expression) %>%
    left_join(drug_gene_names, by = "EnsGeneID") %>%
    relocate(HugoID)     
}
```

``` r
# calculating median expression of druggable genes per disease
aRMS_polyA_drug_medians <- compute_gene_medians_polyA(aRMS_polyA_drug) %>%
  mutate(Disease = "aRMS", Compendia = "PolyA", name = "polyA_medians")

aRMS_riboD_drug_medians <- compute_gene_medians_riboD(aRMS_riboD_drug) %>%
  mutate(Disease = "aRMS", Compendia = "RiboD", name = "riboD_medians")

SS_polyA_drug_medians <- compute_gene_medians_polyA(SS_polyA_drug) %>%
  mutate(Disease = "SS", Compendia = "PolyA", name = "polyA_medians") 

SS_riboD_drug_medians <- compute_gene_medians_riboD(SS_riboD_drug) %>%
  mutate(Disease = "SS", Compendia = "RiboD", name = "riboD_medians") 

WT_polyA_drug_medians <- compute_gene_medians_polyA(WT_polyA_drug) %>%
  mutate(Disease = "WT", Compendia = "PolyA", name = "polyA_medians") 

WT_riboD_drug_medians <- compute_gene_medians_riboD(WT_riboD_drug) %>%
  mutate(Disease = "WT", Compendia = "RiboD", name = "riboD_medians")

NB_polyA_drug_medians <- compute_gene_medians_polyA(NB_polyA_drug) %>%
  mutate(Disease = "NB", Compendia = "PolyA", name = "polyA_medians") 

NB_riboD_drug_medians <- compute_gene_medians_riboD(NB_riboD_drug) %>%
  mutate(Disease = "NB", Compendia = "RiboD", name = "riboD_medians") 

ALL_polyA_drug_medians <- compute_gene_medians_polyA(ALL_polyA_drug) %>%
  mutate(Disease = "ALL", Compendia = "PolyA", name = "polyA_medians") 

ALL_riboD_drug_medians <- compute_gene_medians_riboD(ALL_riboD_drug) %>%
  mutate(Disease = "ALL", Compendia = "RiboD", name = "riboD_medians") 

AML_polyA_drug_medians <- compute_gene_medians_polyA(AML_polyA_drug) %>%
  mutate(Disease = "AML", Compendia = "PolyA", name = "polyA_medians") 

AML_riboD_drug_medians <- compute_gene_medians_riboD(AML_riboD_drug) %>%
  mutate(Disease = "AML", Compendia = "RiboD", name = "riboD_medians") 
```

``` r
combined_drug_medians <- bind_rows(
  aRMS_polyA_drug_medians,
  aRMS_riboD_drug_medians,
  SS_polyA_drug_medians,
  SS_riboD_drug_medians,
  NB_polyA_drug_medians,
  NB_riboD_drug_medians,
  WT_polyA_drug_medians,
  WT_riboD_drug_medians,
  ALL_polyA_drug_medians,
  ALL_riboD_drug_medians,
  AML_polyA_drug_medians,
  AML_riboD_drug_medians
)
# write_tsv(combined_drug_medians, "../../output_data/Fig1H/combined_drug_medians.tsv.gz")
```

``` r
# Define consistent color-blind–safe palette
scale_color_compendia <- function() {
  scale_color_manual(values = c(
    "RiboD" = "#E69F00",  # yellow
    "PolyA" = "#0072B2"   # blue
  ))
}
```

``` r
# Loop through diseases and generate plots
Fig1H_SS_ensembl <- map(unique(combined_drug_medians$Disease), function(d) {

 df <- combined_drug_medians %>% 
  filter(Disease == d)

# order ensemblIDs by median expression
df <- combined_drug_medians %>% 
  filter(Disease == d) %>%
  mutate(EnsGeneID = factor(EnsGeneID, levels = df %>% 
                         group_by(EnsGeneID) %>% 
                         summarize(median_expr = median(Expression), .groups = "drop") %>% 
                         arrange(median_expr) %>% 
                         pull(EnsGeneID)))

  ggplot(df, aes(x = EnsGeneID, y = Expression, color = Compendia)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # 1. vertical lines for each gene 
    geom_vline(aes(xintercept = as.numeric(EnsGeneID)), color = "grey85", linewidth = 0.3) + 
    # 2. points that stay centered on those lines 
    geom_point(size = 1.4, alpha = 0.75) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      title = paste(d),
      x = "Treehouse Druggable Genes",
      # y = "Expression log2(TPM+1)",
      color = "Library Prep"
    ) +
    ylab(bquote(Median~log[2](TPM+1))) +
    theme(
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(), 
      # axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 12),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 15), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "bottom"
    )
})

# Assign names for easy access
names(Fig1H_SS_ensembl) <- unique(combined_drug_medians$Disease)

# View plots
Fig1H_SS_ensembl$SS
```

![](Fig1H_files/figure-commonmark/Fig1H_SS-1.png)

``` r
ggsave("../../Figures/Fig1H_SS_ensembl.png", Fig1H_SS_ensembl$SS, width = 11)
```

    Saving 11 x 5 in image

``` r
# ggsave("../../Figures/Fig1H.tif", Fig1H_SS$SS, width = 11)
```

Fig S7

``` r
# Loop through diseases and generate plots
Fig1H_ensembl <- map(unique(combined_drug_medians$Disease), function(d) {

 df <- combined_drug_medians %>%
  filter(Disease == d)

df <- combined_drug_medians %>%
  filter(Disease == d) %>%
  mutate(EnsGeneID = factor(EnsGeneID, levels = df %>%
                         group_by(EnsGeneID) %>%
                         summarize(median_expr = median(Expression), .groups = "drop") %>%
                         arrange(median_expr) %>%
                         pull(EnsGeneID)))

  ggplot(df, aes(x = EnsGeneID, y = Expression, color = Compendia)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # 1. vertical lines for each gene
    geom_vline(aes(xintercept = as.numeric(EnsGeneID)), color = "grey85", linewidth = 0.3) +
    # 2. points that stay centered on those lines
    geom_point(size = 3, alpha = 0.75) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      title = paste(d),
      x = NULL,
      # y = "Expression log2(TPM+1)",
      color = "Library Prep"
    ) +
    ylab(bquote(Median~log[2](TPM+1))) +
    theme(
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(),
      # axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 28),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 30, face = "bold"),
      axis.title.x = element_text(angle = 0, hjust = 1, size = 30),
      strip.text = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
      legend.position = "bottom",
      legend.text = element_text(face = "bold", size = 26),
      legend.title = element_text(face = "bold", size = 26)
    )
})

# Assign names for easy access
names(Fig1H_ensembl) <- unique(combined_drug_medians$Disease)

# View plots
Fig1H_ensembl
```

    $aRMS

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-1.png)


    $SS

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-2.png)


    $NB

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-3.png)


    $WT

    Warning: Removed 1 row containing missing values or values outside the scale range
    (`geom_point()`).

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-4.png)


    $ALL

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-5.png)


    $AML

![](Fig1H_files/figure-commonmark/Fig1H_ensembl-6.png)

``` r
FigS7_ensembl <- wrap_plots(Fig1H_ensembl, ncol = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") +
  plot_annotation(
    theme = theme(legend.position = "bottom")
  ) &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 5, r = 5)
    )

FigS7_ensembl
```

    Warning: Removed 1 row containing missing values or values outside the scale range
    (`geom_point()`).

![](Fig1H_files/figure-commonmark/FigS7_ensembl-1.png)

``` r
# ggsave("../../Figures/FigS7_ensembl.png", FigS7_ensembl, width = 30, height = 40, dpi = 300)
# ggsave("../../Figures/FigS7_ensembl.png.tif", FigS7_ensembl.png, width = 36, height = 36, dpi = 300, scale = 0.75)
```

``` r
# labeling x-axis based on hugoID
# Loop through diseases and generate plots
Fig1H_SS_hugo <- map(unique(combined_drug_medians$Disease), function(d) {

 df <- combined_drug_medians %>% 
  filter(Disease == d)

# order ensemblIDs by median expression
  gene_order <- df %>% 
    group_by(EnsGeneID) %>%
    summarize(median_expr = median(Expression), .groups = "drop") %>%
    arrange(median_expr) %>%
    pull(EnsGeneID)

  df <- df %>% mutate(EnsGeneID = factor(EnsGeneID, levels = gene_order))
  
  # named lookup: names = EnsGeneID, values = HugoID
  label_lookup <- df %>%
    distinct(EnsGeneID, HugoID) %>%
    { setNames(.$HugoID, as.character(.$EnsGeneID)) }

  ggplot(df, aes(x = EnsGeneID, y = Expression, color = Compendia)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # 1. vertical lines for each gene 
    geom_vline(aes(xintercept = as.numeric(EnsGeneID)), color = "grey85", linewidth = 0.3) + 
    # 2. points that stay centered on those lines 
    geom_point(size = 1.4, alpha = 0.75) +
    scale_x_discrete(
      labels = label_lookup,
      expand = expansion(mult = c(0.01, 0.01))
      ) +
    scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      title = paste(d),
      x = "Treehouse Druggable Genes",
      # y = "Expression log2(TPM+1)",
      color = "Library Prep"
    ) +
    ylab(bquote(Median~log[2](TPM+1))) +
    theme(
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(), 
      # axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 12),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 15), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "bottom"
    )
})

# Assign names for easy access
names(Fig1H_SS_hugo) <- unique(combined_drug_medians$Disease)

# View plots
Fig1H_SS_hugo$SS
```

![](Fig1H_files/figure-commonmark/Fig1H_SS_hugo-1.png)

``` r
# ggsave("../../Figures/Fig1H_SS_hugo.png", Fig1H_SS_hugo$SS, width = 11)
```

``` r
# labeling x-axis based on hugoID
# Loop through diseases and generate plots
Fig1H_hugo <- map(unique(combined_drug_medians$Disease), function(d) {

 df <- combined_drug_medians %>% 
  filter(Disease == d)

# order ensemblIDs by median expression
  gene_order <- df %>% 
    group_by(EnsGeneID) %>%
    summarize(median_expr = median(Expression), .groups = "drop") %>%
    arrange(median_expr) %>%
    pull(EnsGeneID)

  df <- df %>% mutate(EnsGeneID = factor(EnsGeneID, levels = gene_order))
  
  # named lookup: names = EnsGeneID, values = HugoID
  label_lookup <- df %>%
    distinct(EnsGeneID, HugoID) %>%
    { setNames(.$HugoID, as.character(.$EnsGeneID)) }

  ggplot(df, aes(x = EnsGeneID, y = Expression, color = Compendia)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # 1. vertical lines for each gene 
    geom_vline(aes(xintercept = as.numeric(EnsGeneID)), color = "grey85", linewidth = 0.3) + 
    # 2. points that stay centered on those lines 
    geom_point(size = 3, alpha = 0.75) +
    scale_x_discrete(
      labels = label_lookup,
      expand = expansion(mult = c(0.01, 0.01))
      ) +
    scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      title = paste(d),
      x = NULL,
      # y = "Expression log2(TPM+1)",
      color = "Library Prep"
    ) +
    ylab(bquote(Median~log[2](TPM+1))) +
    theme(
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(), 
      # axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 28),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 28),
      axis.title.y = element_text(angle = 90, hjust = 0.5, size = 30, face = "bold"),
      axis.title.x = element_text(angle = 0, hjust = 1, size = 30),
      strip.text = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
      legend.position = "bottom",
      legend.text = element_text(face = "bold", size = 26),
      legend.title = element_text(face = "bold", size = 26)
    )
})

# Assign names for easy access
names(Fig1H_hugo) <- unique(combined_drug_medians$Disease)

# View plots
Fig1H_hugo
```

    $aRMS

![](Fig1H_files/figure-commonmark/Fig1H_hugo-1.png)


    $SS

![](Fig1H_files/figure-commonmark/Fig1H_hugo-2.png)


    $NB

![](Fig1H_files/figure-commonmark/Fig1H_hugo-3.png)


    $WT

    Warning: Removed 1 row containing missing values or values outside the scale range
    (`geom_point()`).

![](Fig1H_files/figure-commonmark/Fig1H_hugo-4.png)


    $ALL

![](Fig1H_files/figure-commonmark/Fig1H_hugo-5.png)


    $AML

![](Fig1H_files/figure-commonmark/Fig1H_hugo-6.png)

``` r
FigS7_hugo <- wrap_plots(Fig1H_hugo, ncol = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") +
  plot_annotation(
    theme = theme(legend.position = "bottom")
  ) &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 5, r = 5)
    )

FigS7_hugo
```

    Warning: Removed 1 row containing missing values or values outside the scale range
    (`geom_point()`).

![](Fig1H_files/figure-commonmark/FigS7_hugo-1.png)

``` r
# ggsave("../../Figures/FigS7_hugo.png", FigS7_hugo, width = 30, height = 40, dpi = 300)
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
     date     2026-10-06
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     package      * version date (UTC) lib source
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
     bit64          4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
     cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
     cowplot      * 1.2.0   2025-07-07 [1] CRAN (R 4.5.0)
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
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
     knitr          1.50    2025-03-16 [1] CRAN (R 4.5.0)
     labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
     lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
     lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
     magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
     patchwork    * 1.3.2   2025-08-25 [1] CRAN (R 4.5.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
     purrr        * 1.1.0   2025-07-10 [1] CRAN (R 4.5.0)
     R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
     ragg           1.5.0   2025-09-02 [1] CRAN (R 4.5.0)
     RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
     readr        * 2.1.5   2024-01-10 [1] CRAN (R 4.5.0)
     rlang          1.2.0   2026-04-06 [1] CRAN (R 4.5.2)
     rmarkdown      2.30    2025-09-28 [1] CRAN (R 4.5.0)
     rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
     S7             0.2.0   2024-11-07 [1] CRAN (R 4.5.0)
     scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
     sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
     stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
     stringr      * 1.5.2   2025-09-08 [1] CRAN (R 4.5.0)
     systemfonts    1.3.1   2025-10-01 [1] CRAN (R 4.5.0)
     textshaping    1.0.4   2025-10-10 [1] CRAN (R 4.5.0)
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
