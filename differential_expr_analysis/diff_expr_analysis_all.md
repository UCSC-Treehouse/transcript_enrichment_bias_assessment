# diff_expr_analysis_all


## **Fig3**

**Representative assessment of biases introduced in differential
expression analysis when comparing two diseases across transcript
enrichment methods**

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
library(ggplot2)
library(cowplot)
```


    Attaching package: 'cowplot'

    The following object is masked from 'package:lubridate':

        stamp

``` r
library(VennDiagram)
```

    Loading required package: grid
    Loading required package: futile.logger

``` r
aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../output_data/aRMS_SS/aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased.rds") #for fig 2A - venn diagram

aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../output_data/aRMS_SS/aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased.rds") #for fig 2B - venn diagram

aRMS_SS_up_polyAbiased <- read_rds("../output_data/aRMS_SS/aRMS_SS_up_polyAbiased.rds") #for fig 2C - histogram

aRMS_SS_up_riboDbiased <- read_rds("../output_data/aRMS_SS/aRMS_SS_up_riboDbiased.rds") #for fig 2D - histogram

up_polyAbiased_all <- read_rds("../output_data/all_disease_comparisons/up_polyAbiased_all.rds") #for fig 2E - venn diagram

up_riboDbiased_all <- read_rds("../output_data/all_disease_comparisons/up_riboDbiased_all.rds") #for fig 2E - venn diagram
```

``` r
# for fig 3A
Fig_3A <- venn.diagram(
  x = aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased,
  category.names = c(
    "PolyA unbiased",
    "PolyA biased",
    "RiboD biased"
    ),
  filename = NULL, # Save as a PNG file
  output = TRUE,
  print.mode = c("raw", "percent"),
  # Customize appearance (optional)
  fill = c("#BB5566", "#0072B2", "#E69F00"),
#  cat.col = c("#BB5566", "#0072B2", "#E69F00"),
  cex = 1, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  # main = "Upregulated genes in SS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-10-27 16:03:41] $x
    INFO [2026-10-27 16:03:41] aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $category.names
    INFO [2026-10-27 16:03:41] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $filename
    INFO [2026-10-27 16:03:41] NULL
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $output
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $print.mode
    INFO [2026-10-27 16:03:41] c("raw", "percent")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fill
    INFO [2026-10-27 16:03:41] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cex
    INFO [2026-10-27 16:03:41] [1] 1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.cex
    INFO [2026-10-27 16:03:41] [1] 0.9
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.dist
    INFO [2026-10-27 16:03:41] c(0.04, 0.04, 0.04)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.default.pos
    INFO [2026-10-27 16:03:41] [1] "outer"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.pos
    INFO [2026-10-27 16:03:41] c(-14, 14, 175)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.cex
    INFO [2026-10-27 16:03:41] [1] 1.1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $disable.logging
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $sigdigs
    INFO [2026-10-27 16:03:41] [1] 3
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 

``` r
Fig_3A
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3A-1.png)

``` r
# for fig 3B
Fig_3B <- venn.diagram(
  x = aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased,
  category.names = c(
    "RiboD unbiased",
    "PolyA biased",
    "RiboD biased"
    ),
  filename = NULL, # Save as a PNG file
  output = TRUE,
  print.mode = c("raw", "percent"),
  # Customize appearance (optional)
  fill = c("#BB5566", "#0072B2", "#E69F00"),
#  cat.col = c("#BB5566", "#0072B2", "#E69F00"),
  cex = 1, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 2000,
  # width = 2000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  # main = "Upregulated genes in SS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-10-27 16:03:41] $x
    INFO [2026-10-27 16:03:41] aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $category.names
    INFO [2026-10-27 16:03:41] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $filename
    INFO [2026-10-27 16:03:41] NULL
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $output
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $print.mode
    INFO [2026-10-27 16:03:41] c("raw", "percent")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fill
    INFO [2026-10-27 16:03:41] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cex
    INFO [2026-10-27 16:03:41] [1] 1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.cex
    INFO [2026-10-27 16:03:41] [1] 0.9
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.dist
    INFO [2026-10-27 16:03:41] c(0.04, 0.04, 0.04)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.default.pos
    INFO [2026-10-27 16:03:41] [1] "outer"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.pos
    INFO [2026-10-27 16:03:41] c(-14, 14, 175)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.cex
    INFO [2026-10-27 16:03:41] [1] 1.1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $disable.logging
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $sigdigs
    INFO [2026-10-27 16:03:41] [1] 3
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 

``` r
Fig_3B
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3B-1.png)

``` r
aRMS_SS_up_polyAbiased_median <- aRMS_SS_up_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
aRMS_SS_up_polyAbiased_median
```

| median_value |
|-------------:|
|     1.322192 |

``` r
# for fig 3C
Fig_3C <- ggplot(aRMS_SS_up_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
#    title = "Unique to PolyA Biased",
#    subtitle = "polyA biased, Median LFC = 1.32",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = aRMS_SS_up_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=3.7, y=580, label= "Median = 1.32")

Fig_3C
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3C-1.png)

``` r
aRMS_SS_up_riboDbiased_median <- aRMS_SS_up_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
aRMS_SS_up_riboDbiased_median
```

| median_value |
|-------------:|
|     2.761514 |

``` r
# for fig 3D
Fig_3D <- ggplot(aRMS_SS_up_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
#    title = "Unique to RiboD Biased",
    # title = "Distribution of log2FoldChange values of \n genes uniquely upregulated in SS riboD relative to aRMS polyA",
    # subtitle = "riboD biased, Median LFC = 2.76",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = aRMS_SS_up_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=7, y=1200, label= "Median = 2.76")

Fig_3D 
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3D-1.png)

``` r
# for fig3E
Fig_3E <- venn.diagram(
  x = up_polyAbiased_all,
  category.names = c(
    "SS to aRMS",
    "WT to NB",
    "AML to ALL"
    ),
  filename = NULL, # Save as a rds file
  output = TRUE,
  print.mode = c("raw", "percent"),
  sigdigs = 3, # the amount of significant digits in percent
  # Customize appearance (optional)
  fill = c("#0072B2", "#0098ed", "#87CEEB"),
#  cat.col = c("#0072B2", "#0072B2", "#0072B2"),
  cex = 1, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 2000,
  # width = 2000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
    cat.fontfamily = "sans",
    fontfamily = "sans",
  main = "PolyA Biased",
  main.fontfamily = "sans",
  main.cex = 1.1, # Font size for main title
  disable.logging = TRUE
)
```

    INFO [2026-10-27 16:03:41] $x
    INFO [2026-10-27 16:03:41] up_polyAbiased_all
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $category.names
    INFO [2026-10-27 16:03:41] c("SS to aRMS", "WT to NB", "AML to ALL")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $filename
    INFO [2026-10-27 16:03:41] NULL
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $output
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $print.mode
    INFO [2026-10-27 16:03:41] c("raw", "percent")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $sigdigs
    INFO [2026-10-27 16:03:41] [1] 3
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fill
    INFO [2026-10-27 16:03:41] c("#0072B2", "#0098ed", "#87CEEB")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cex
    INFO [2026-10-27 16:03:41] [1] 1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.cex
    INFO [2026-10-27 16:03:41] [1] 0.9
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.dist
    INFO [2026-10-27 16:03:41] c(0.04, 0.04, 0.04)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.default.pos
    INFO [2026-10-27 16:03:41] [1] "outer"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.pos
    INFO [2026-10-27 16:03:41] c(-14, 14, 175)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main
    INFO [2026-10-27 16:03:41] [1] "PolyA Biased"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.cex
    INFO [2026-10-27 16:03:41] [1] 1.1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $disable.logging
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 

``` r
Fig_3E
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3E-1.png)

``` r
# for fig3F
Fig_3F <- venn.diagram(
  x = up_riboDbiased_all,
  category.names = c(
    "SS to aRMS",
    "WT to NB",
    "AML to ALL"
    ),
  filename = NULL, # Save as a rds file
  output = TRUE,
  print.mode = c("raw", "percent"),
  sigdigs = 3, # the amount of significant digits in percent
  # Customize appearance (optional)
  fill = c("#E69F00", "#ff9b00", "#e6b200"),
#  cat.col = c("#E69F00", "#E69F00", "#E69F00"),
  cex = 1.1, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04), # distance of category name to VD
  # height = 2000,
  # width = 2000,
  cat.default.pos = "outer", # location of category names
  cat.pos = c(-14, 14, 175), # position of category names
  main = "RiboD Biased",
    fontfamily = "sans",
  main.fontfamily = "sans",
    cat.fontfamily = "sans",
  main.cex = 1.1, # Font size for main title
  disable.logging = TRUE
)
```

    INFO [2026-10-27 16:03:41] $x
    INFO [2026-10-27 16:03:41] up_riboDbiased_all
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $category.names
    INFO [2026-10-27 16:03:41] c("SS to aRMS", "WT to NB", "AML to ALL")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $filename
    INFO [2026-10-27 16:03:41] NULL
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $output
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $print.mode
    INFO [2026-10-27 16:03:41] c("raw", "percent")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $sigdigs
    INFO [2026-10-27 16:03:41] [1] 3
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fill
    INFO [2026-10-27 16:03:41] c("#E69F00", "#ff9b00", "#e6b200")
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cex
    INFO [2026-10-27 16:03:41] [1] 1.1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.cex
    INFO [2026-10-27 16:03:41] [1] 0.9
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.dist
    INFO [2026-10-27 16:03:41] c(0.04, 0.04, 0.04)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.default.pos
    INFO [2026-10-27 16:03:41] [1] "outer"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.pos
    INFO [2026-10-27 16:03:41] c(-14, 14, 175)
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main
    INFO [2026-10-27 16:03:41] [1] "RiboD Biased"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $cat.fontfamily
    INFO [2026-10-27 16:03:41] [1] "sans"
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $main.cex
    INFO [2026-10-27 16:03:41] [1] 1.1
    INFO [2026-10-27 16:03:41] 
    INFO [2026-10-27 16:03:41] $disable.logging
    INFO [2026-10-27 16:03:41] [1] TRUE
    INFO [2026-10-27 16:03:41] 

``` r
Fig_3F
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3F-1.png)

``` r
Fig_3 <- plot_grid(
  Fig_3A,
  Fig_3B,
  Fig_3C,
  Fig_3D,
  Fig_3E,
  Fig_3F,
  ncol = 2,
  labels = c('A', 'B', 'C', 'D', 'E', 'F'),
  label_size = 15)

ggsave("../Figures/Fig_3.png", Fig_3, width = 8, height = 12, dpi = 300)
# ggsave("../Figures/Fig_3.tif", Fig_3, width = 8, height = 12, dpi = 300)

Fig_3
```

![](diff_expr_analysis_all_files/figure-commonmark/Fig_3-1.png)

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
     date     2026-10-27
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     ! package        * version date (UTC) lib source
     P cli              3.6.5   2025-04-23 [?] CRAN (R 4.5.0)
     P cowplot        * 1.2.0   2025-07-07 [?] RSPM
     P digest           0.6.37  2024-08-19 [?] CRAN (R 4.5.0)
     P dplyr          * 1.2.1   2026-04-03 [?] RSPM
     P evaluate         1.0.5   2025-08-27 [?] RSPM
     P farver           2.1.2   2024-05-13 [?] RSPM
     P fastmap          1.2.0   2024-05-15 [?] RSPM
     P forcats        * 1.0.1   2025-09-25 [?] RSPM
     P formatR          1.14    2023-01-17 [?] CRAN (R 4.5.0)
     P futile.logger  * 1.4.9   2025-12-29 [?] RSPM
     P futile.options   1.0.1   2018-04-20 [?] CRAN (R 4.5.0)
     P generics         0.1.4   2025-05-09 [?] RSPM
     P ggplot2        * 4.0.3   2026-04-22 [?] RSPM
     P glue             1.8.0   2024-09-30 [?] CRAN (R 4.5.0)
     P gtable           0.3.6   2024-10-25 [?] RSPM
     P hms              1.1.4   2025-10-17 [?] RSPM
     P htmltools        0.5.8.1 2024-04-04 [?] CRAN (R 4.5.0)
     P jsonlite         2.0.0   2025-03-27 [?] RSPM
     P knitr            1.50    2025-03-16 [?] CRAN (R 4.5.0)
     P labeling         0.4.3   2023-08-29 [?] RSPM
     P lambda.r         1.2.6   2026-10-13 [?] RSPM
     P lifecycle        1.0.5   2026-01-08 [?] RSPM
     P lubridate      * 1.9.5   2026-02-04 [?] RSPM
     P magrittr         2.0.5   2026-04-04 [?] RSPM
     P pillar           1.11.1  2025-09-17 [?] RSPM
     P pkgconfig        2.0.3   2019-09-22 [?] RSPM
     P purrr          * 1.2.2   2026-04-10 [?] RSPM
     P R6               2.6.1   2025-02-15 [?] RSPM
     P ragg             1.5.2   2026-03-23 [?] RSPM
     P RColorBrewer     1.1-3   2022-04-03 [?] RSPM
     P readr          * 2.2.0   2026-02-19 [?] RSPM
       renv             1.1.5   2025-07-24 [1] CRAN (R 4.5.0)
     P rlang            1.2.0   2026-04-06 [?] RSPM
     P rmarkdown        2.30    2025-09-28 [?] CRAN (R 4.5.0)
     P rstudioapi       0.19.0  2026-06-11 [?] RSPM
     P S7               0.2.2   2026-04-22 [?] RSPM
     P scales           1.4.0   2025-04-24 [?] RSPM
     P sessioninfo      1.2.4   2026-06-04 [?] RSPM
     P stringi          1.8.9   2026-08-04 [?] RSPM
     P stringr        * 1.6.0   2025-11-04 [?] RSPM
     P systemfonts      1.3.2   2026-03-05 [?] RSPM
     P textshaping      1.0.5   2026-03-06 [?] RSPM
     P tibble         * 3.3.1   2026-01-11 [?] RSPM
     P tidyr          * 1.3.2   2025-12-19 [?] RSPM
     P tidyselect       1.2.1   2024-03-11 [?] RSPM
     P tidyverse      * 2.0.0   2023-02-22 [?] RSPM
     P timechange       0.4.0   2026-01-29 [?] RSPM
     P tzdb             0.5.0   2025-03-15 [?] RSPM
     P vctrs            0.7.3   2026-04-11 [?] RSPM
     P VennDiagram    * 1.8.2   2026-01-11 [?] RSPM
     P withr            3.0.3   2026-06-19 [?] RSPM
     P xfun             0.55    2025-12-16 [?] CRAN (R 4.5.2)
     P yaml             2.3.10  2024-07-26 [?] CRAN (R 4.5.0)

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/differential_expr_analysis/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Users/maryke/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
