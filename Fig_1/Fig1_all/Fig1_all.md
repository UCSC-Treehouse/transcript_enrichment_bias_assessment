# Fig1_all


## Fig1

**Gene and sample variability from RNA-seq libraries prepared from
different transcript enrichment methods.**

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
library(ggplot2)
library(patchwork)
library(magick)
```

    Linking to ImageMagick 6.9.13.29
    Enabled features: cairo, fontconfig, freetype, heic, lcms, pango, raw, rsvg, webp
    Disabled features: fftw, ghostscript, x11

``` r
library(rsvg)
```

    Linking to librsvg 2.61.0

``` r
library(ggplotify)
library(grid)
library(ggrepel) 
library(cowplot)
```


    Attaching package: 'cowplot'

    The following object is masked from 'package:patchwork':

        align_plots

    The following object is masked from 'package:lubridate':

        stamp

``` r
Fig1A <- image_read("../../Figures/UMAPs/Compendia_Type.svg") %>% 
  image_ggplot()
Fig1B <- image_read("../../Figures/UMAPs/10_Most_Sampled_Diseases.svg") %>% 
  image_ggplot()
```

### Fig1C

``` r
all_medians_bins_1C <- read_tsv("../../output_data/Fig1C/all_medians_bins.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 725976 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): Gene, Disease, Compendia, ExprBin
    dbl (1): median_expr

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
bin_prop_1C <- read_tsv("../../output_data/Fig1C/bin_prop.tsv.gz") %>%
  mutate(ExprBin = factor(ExprBin, levels = c("0", "0-0.09", "0.1-0.99", "1-2.99", "3-4.99", ">5"))) %>%
  filter(Disease == "SS")
```

    Rows: 72 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (3): Compendia, ExprBin, Disease
    dbl (2): n, Proportion

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# expr_bins <- c(0, 1e-10, 0.09, 0.99, 2.99, 4.99, Inf)
# expr_labels <- c("0","0-0.09", "0.1-0.99", "1-2.99", "3-4.99", ">5")
# bin_prop_1C <- bin_prop_1C %>%
#   mutate(ExprBin = factor(ExprBin, levels = c("0", "0-0.09", "0.1-0.99", "1-2.99", "3-4.99", ">5"))) %>%
#   filter(Disease == "SS")
```

``` r
theme_1C <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "right",
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
    "polyA" = "#0072B2",  # blue
    "riboD" = "#E69F00"   # yellow
  ))
}
scale_color_compendia <- function() {
  scale_color_manual(values = c( 
    "polyA" = "#0072B2",  # blue
    "riboD" = "#E69F00"   # yellow
    ))
  }
```

``` r
Fig1C_SS <- ggplot(bin_prop_1C, aes(x = Compendia, y = Proportion, fill = ExprBin)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis_d(option = "D", name = "Expression \nin log2(TPM+1)") +
    labs(
      # x = "Library Prep",
      y = "Proportion of Genes",
      fill = "Expression Level"
    ) +
    theme_1C()

Fig1C_SS
```

![](Fig1_all_files/figure-commonmark/Fig1C-1.png)

### Fig1D

``` r
sig_bins_all_1D <- read_tsv("../../output_data/Fig1D/sig_bins_all.tsv.gz") %>%
  filter(Disease == "SS")  %>%
    select(Bin, stars_adj)
```

    Rows: 36 Columns: 8
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): Disease, Bin, test_used, stars_adj
    dbl (2): p_value, padj
    lgl (2): polyA, riboD

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_counts_1D <- read_tsv("../../output_data/Fig1D/combined_counts.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 1200 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): Sample, Bin, Disease, Compendia
    dbl (1): Count

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all_stats_1D <- read_tsv("../../output_data/Fig1D/all_stats.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 36 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): Disease, Bin
    dbl (4): Shapiro_polyA_p, Shapiro_riboD_p, VarTest_p, TTest_p

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
theme_1d <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
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
```

``` r
facet_labels <- c(
  "0" = "x = 0",
  "0-0.09" = "0 > x > 0.1",
  "0.1-0.99" = "0.1 > x > 1",
  "1-2.99" = "1 > x > 3",
  "3-4.99" = "3 > x > 5",
  ">5" = "x > 5"
)

# List of diseases
diseases <- unique(combined_counts_1D$Disease)

max_y <- max(combined_counts_1D$Count, na.rm = TRUE)
```

``` r
# y-position for stars
y_star <- max(combined_counts_1D$Count, na.rm = TRUE) * 1.05

Fig1D_SS <- ggplot(combined_counts_1D, aes(x = Compendia, y = Count, fill = Compendia)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6,
                 position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = Compendia),
                position = position_jitter(width = 0.15),
                size = 1.5, alpha = 0.7) +

    # ---- ADD STARS ABOVE EACH FACET ----
    geom_text(
      data = sig_bins_all_1D,
      aes(x = 1.5, y = y_star, label = stars_adj),   # x = 1.5 centers between polyA/riboD
      inherit.aes = FALSE,
      size = 6
    ) +

    facet_grid(. ~ Bin, labeller = labeller(Bin = facet_labels)) +
    scale_fill_compendia() +
    scale_color_compendia() +
#    coord_cartesian(ylim = c(0, y_star * 1.1)) +
    labs(
      x = "log2(TPM+1) = x",
      y = "Number of Genes Expressed"
    ) +
    theme_minimal() +
    theme_1d() +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(1, "lines")
    )

Fig1D_SS
```

![](Fig1_all_files/figure-commonmark/Fig1D-1.png)

### Fig1E

``` r
median_gene_map <- read_tsv("../../output_data/Fig1E_F/median_gene_map_df.tsv.gz")
```

    Rows: 1 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (6): SS, aRMS, WT, NB, ALL, AML

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_long <- read_tsv("../../output_data/Fig1E_F/combined_long.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 400 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): Gene, TH_ID, Disease, Compendia
    dbl (1): Expression

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sig_results <- read_tsv("../../output_data/Fig1E_F/sig_results.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 12 Columns: 9
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (5): Disease, Gene, test_used, stars, stars_adj
    dbl (2): p_value, padj
    lgl (2): polyA, riboD

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
color_theme_1E <- function(base_size = 14) {
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
```

``` r
median_gene <- "MTMR2"


  # SIGNIFICANCE STARS
star_df <- tibble(
    Disease = "SS",
    Gene = c("HIST1H1B", median_gene),
    PlotType = c("HIST1H1B (PolyA-)", "Median Ratio Gene (PolyA+)"),
    stars_adj = sig_results %>%
      filter(Disease == "SS", Gene %in% c("HIST1H1B", median_gene)) %>%
      arrange(match(Gene, c("HIST1H1B", median_gene))) %>%
      pull(stars_adj)
  )

# y-position for stars (slightly above max)
y_star <- max(combined_long$Expression, na.rm = TRUE) * 1.05

Fig1E_SS <- ggplot(combined_long, aes(x = Gene, y = Expression, fill = Compendia)) +
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
      size = 12
    ) +
    scale_fill_compendia() +
    scale_color_compendia() +
    coord_cartesian(ylim = c(0, 10)) +
    # facet_wrap(~ PlotType, ncol = 2, scales = "free_x") +
    labs(
      x = NULL,
      y = "Expression log2(TPM+1)"
    ) +
    color_theme_1E() +
    theme(
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 20),
      strip.text = element_blank(),
      plot.title = element_text(vjust = -2),
      legend.position = "bottom",
      legend.text = element_text(size = 16),       # <-- label text size
      legend.title = element_text(size = 18),      # <-- title text size
      legend.key.size = unit(1, "cm"), # <-- box/key size
      legend.margin = margin(t = -2) 
    ) 
Fig1E_SS
```

![](Fig1_all_files/figure-commonmark/Fig1E_SS-1.png)

### Fig1F

``` r
ss_ratios <- read_tsv("../../output_data/Fig1E_F/ss_ratios.tsv.gz")
```

    Rows: 60498 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): Gene, Disease
    dbl (1): Ratio

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SS_most_average_gene_hugo <- read_tsv("../../output_data/Fig1E_F/SS_most_average_gene_hugo.tsv.gz")
```

    Rows: 1 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (1): Gene
    dbl (4): SS_polyA_median, SS_riboD_median, SS_median_ratio, dist

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SS_HIST1H1B <- read_tsv("../../output_data/Fig1E_F/SS_HIST1H1B.tsv.gz")
```

    Rows: 1 Columns: 4
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (1): Gene
    dbl (3): SS_polyA_median, SS_riboD_median, SS_median_ratio

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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
# Fig1F - SS standalone
Fig1F_SS <- ggplot(ss_ratios, aes(x = Ratio)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  
  geom_vline(xintercept = SS_most_average_gene_hugo$SS_median_ratio, color = "#0072B2", linetype = "dashed", linewidth = 1) +
  
  geom_text(aes(x = 2.5, y = 2000, label = SS_most_average_gene_hugo$Gene), angle = 0, vjust = -0.4, hjust = 0.3, color = "#0072B2", size = 6) +
  
  geom_vline(xintercept = SS_HIST1H1B$SS_median_ratio, color = "#E69F00", linetype = "dashed", linewidth = 1) +
  
  geom_text(aes(x = 0.9, y = 2000, label = SS_HIST1H1B$Gene), angle = 0, vjust = -5, hjust = 0.2, color = "#E69F00", size = 6) +
  
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 4000)) +
  labs(
    title = "SS",
    x = "Median Expression Ratio",
    y = "Number of Genes"
  ) +
  theme_Fig1F() +
      theme(
      axis.title.y = element_text(size = 22),
      axis.title.x = element_text(size = 22),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      strip.text = element_blank(),
      plot.title = element_text(vjust = -2)
    )

Fig1F_SS
```

    Warning in geom_text(aes(x = 2.5, y = 2000, label = SS_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.9, y = 2000, label = SS_HIST1H1B$Gene), angle = 0, : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 29257 rows containing non-finite outside the scale range
    (`stat_bin()`).

![](Fig1_all_files/figure-commonmark/Fig1F-1.png)

### Fig1G

``` r
combined_medians_hugo_res <- read_tsv("../../output_data/Fig1G/combined_medians_hugo_res.tsv.gz") %>%
  filter(Disease == "SS")
```

    Warning: One or more parsing issues, call `problems()` on your data frame for details,
    e.g.:
      dat <- vroom(...)
      problems(dat)

    Rows: 362988 Columns: 7
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (3): Gene, Disease, label_color
    dbl (3): polyA_medians, riboD_medians, residual
    lgl (1): label

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
poly_gene <- "MTMR2"
nonpoly_gene <- "HIST1H1B"

combined_medians_hugo_res <- combined_medians_hugo_res %>%
  mutate(
    label = case_when(
      Gene == poly_gene[Disease] ~ Gene,
      Gene == nonpoly_gene ~ Gene,            # HIST1H1B + median ratio gene
      TRUE ~ NA_character_
    )
  )
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

diseases <- unique(combined_medians_hugo_res$Disease)
```

``` r
Fig1G_SS <- ggplot(combined_medians_hugo_res, aes(x = polyA_medians, y = riboD_medians)) +
    
   # background points
    geom_point(data = combined_medians_hugo_res %>% filter(is.na(label)),
               color = "grey70", alpha = 0.5, size = 1.2) +

    # highlighted labeled points
    geom_point(data = combined_medians_hugo_res %>% filter(!is.na(label)),
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
      x = "Average PolyA Expression (log2(TPM+1))",
      y = "Average RiboD Expression (log2(TPM+1))",
    ) +
    theme_minimal(base_size = 14) +
    plot_theme()
Fig1G_SS
```

    Warning: Removed 60497 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

![](Fig1_all_files/figure-commonmark/Fig1G-1.png)

### Fig1H

``` r
combined_medians_drug <- read_tsv("../../output_data/Fig1H/combined_medians_drug.tsv.gz") %>%
  filter(Disease == "SS")
```

    Rows: 1356 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (4): Gene, name, Disease, Compendia
    dbl (1): Expression

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_medians_drug <- combined_medians_drug %>%
  mutate(Gene = factor(Gene, levels = combined_medians_drug %>% 
                         group_by(Gene) %>% 
                         summarize(median_expr = median(Expression), .groups = "drop") %>% 
                         arrange(median_expr) %>% 
                         pull(Gene)))
```

``` r
scale_color_compendia <- function() {
  scale_color_manual(values = c(
    "RiboD" = "#E69F00",  # yellow
    "PolyA" = "#0072B2"   # blue
  ))
}

Fig1H_SS <- ggplot(combined_medians_drug, aes(x = Gene, y = Expression, color = Compendia)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # 1. vertical lines for each gene 
    geom_vline(aes(xintercept = as.numeric(Gene)), color = "grey85", linewidth = 0.3) + 
    # 2. points that stay centered on those lines 
    geom_point(size = 1.4, alpha = 0.75) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 11), expand = expansion(mult = c(0, 0.05))) +
    scale_color_compendia() +
    labs(
      x = "Treehouse Druggable Genes",
      y = "Expression log2(TPM+1)",
      color = "Library Prep"
    ) +
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
Fig1H_SS
```

![](Fig1_all_files/figure-commonmark/Fig1H-1.png)

### Fig1 merged

``` r
Fig1_all <- plot_grid(
  NULL,
  NULL,
  Fig1C_SS,
  Fig1D_SS,
  Fig1E_SS,
  Fig1F_SS,
  Fig1G_SS,
  Fig1H_SS,
  ncol = 2,
  label_size = 20,
  labels = "AUTO", 
  label_x = 0.09,     # Nudges slightly away from the absolute left border
  # label_y = 0.95,     # Nudges slightly down from the absolute top border
  hjust = 0,          # Locks alignment to the left edge of the text
  vjust = 1
)
```

    Warning in geom_text(aes(x = 2.5, y = 2000, label = SS_most_average_gene_hugo$Gene), : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning in geom_text(aes(x = 0.9, y = 2000, label = SS_HIST1H1B$Gene), angle = 0, : All aesthetics have length 1, but the data has 60498 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

    Warning: Removed 29257 rows containing non-finite outside the scale range
    (`stat_bin()`).

    Warning: Removed 60497 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).

``` r
Fig1_all
```

![](Fig1_all_files/figure-commonmark/Fig1_all-1.png)

``` r
ggsave("../../Figures/Fig1_all.png", Fig1_all, width = 15, height = 15, dpi = 300)
```
