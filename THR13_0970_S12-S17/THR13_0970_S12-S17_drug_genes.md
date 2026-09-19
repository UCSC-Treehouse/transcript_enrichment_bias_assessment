# THR13_0970_S12-S17_drug_genes


## THR13_0970 matched PolyA / RiboD

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
rsem_log2TPM1_THR13 <- read_tsv("../input_data/matched_THR13_0970/rsem_ensembl_log2TPM1_THR13_0970_S12-S17.tsv.gz")
```

    Rows: 362988 Columns: 8
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (5): full_sample, ensembl_gene_ID, HugoID, lib_prep, sample
    dbl (2): TPM, log2TPM1
    lgl (1): is_treehouse_druggable_gene

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# treehouse druggable gene list
drug_genes <- read_tsv("../input_data/treehouseDruggableGenes_2019-06-12.txt") %>%
  select(gene)
```

    Rows: 115 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): gene, group

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# separating PolyA and RiboD
THR13_polyA <- rsem_log2TPM1_THR13 %>%
  filter(lib_prep == "polyA") %>%
  mutate(Disease = "DIPG4")

THR13_polyA_wide <- THR13_polyA %>%
  select(full_sample, ensembl_gene_ID, HugoID, log2TPM1) %>%
  pivot_wider(names_from = full_sample, values_from = log2TPM1)

THR13_riboD <- rsem_log2TPM1_THR13 %>%
  filter(lib_prep == "riboD") %>%
  mutate(Disease = "DIPG4")

THR13_riboD_wide <- THR13_riboD %>%
  select(full_sample, ensembl_gene_ID, HugoID, log2TPM1) %>%
  pivot_wider(names_from = full_sample, values_from = log2TPM1)
```

``` r
# Find gene medians in polyA
THR13_polyA_medians <- THR13_polyA_wide %>%
  rowwise() %>%
  mutate(median_expr = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(ensembl_gene_ID, HugoID, median_expr)

# Find gene medians in riboD
THR13_riboD_medians <- THR13_riboD_wide %>%
  rowwise() %>%
  mutate(median_expr = median(c_across(starts_with("T")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(ensembl_gene_ID, HugoID, median_expr)
```

``` r
THR13_polyA_medians_drug <- THR13_polyA_medians %>%
  filter(HugoID %in% drug_genes$gene) %>%
  mutate(Disease = "DIPGIV", lib_prep = "PolyA")

THR13_riboD_medians_drug <- THR13_riboD_medians %>%
  filter(HugoID %in% drug_genes$gene) %>%
  mutate(Disease = "DIPGIV", lib_prep = "RiboD")

THR13_drug <- bind_rows(
  THR13_polyA_medians_drug,
  THR13_riboD_medians_drug
)
```

``` r
THR13_drug_form <- THR13_drug %>%
  mutate(HugoID = factor(HugoID, levels = THR13_drug %>%
                           group_by(HugoID) %>%
                           summarize(med_expr = median(median_expr), .groups = "drop") %>%
                           arrange(med_expr) %>%
                           pull(HugoID)))
```

Color theme

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
drug_genes_plot <- ggplot(THR13_drug_form, aes(x = HugoID, y = median_expr, color = lib_prep)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(aes(xintercept = as.numeric(HugoID)), color = "grey85", linewidth = 0.3) + 
  geom_point(size = 1.4, alpha = 0.75) +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      title = "matched DIPGIV",
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
drug_genes_plot
```

![](THR13_0970_S12-S17_drug_genes_files/figure-commonmark/drug_genes_plot-1.png)

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
     date     2026-09-18
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
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
     knitr          1.50    2025-03-16 [1] CRAN (R 4.5.0)
     labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
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
