# Supplemental_Figs


## Fig S8-S15

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
library(patchwork)
```


    Attaching package: 'patchwork'

    The following object is masked from 'package:cowplot':

        align_plots

``` r
# For figS8
aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased.rds")  # venn diagram
WT_NB_up_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_up_polyAunbiased_polyAbiased_riboDbiased.rds") # venn diagram
ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased.rds") # venn diagram

# for figS9
aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased.rds")  # venn diagram
WT_NB_up_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_up_riboDunbiased_polyAbiased_riboDbiased.rds") # venn diagram
ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased.rds") # venn diagram

# for figS10
aRMS_SS_up_polyAbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_up_polyAbiased.rds") # histogram
WT_NB_up_polyAbiased <- read_rds("../../output_data/NB_WT/WT_NB_up_polyAbiased.rds") # histogram
ALL_AML_up_polyAbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_up_polyAbiased.rds") # histogram

# for figS11
aRMS_SS_up_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_up_riboDbiased.rds") # histogram
WT_NB_up_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_up_riboDbiased.rds") # histogram
ALL_AML_up_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_up_riboDbiased.rds") # histogram

# for figS12
aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased.rds")  # venn diagram
WT_NB_down_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_down_polyAunbiased_polyAbiased_riboDbiased.rds") # venn diagram
ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased.rds") # venn diagram

# for figS13
aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased.rds")  # venn diagram
WT_NB_down_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_down_riboDunbiased_polyAbiased_riboDbiased.rds") # venn diagram
ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased.rds") # venn diagram

# for figS14
aRMS_SS_down_polyAbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_down_polyAbiased.rds") # histogram
WT_NB_down_polyAbiased <- read_rds("../../output_data/NB_WT/WT_NB_down_polyAbiased.rds") # histogram
ALL_AML_down_polyAbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_down_polyAbiased.rds") # histogram

# for figS15
aRMS_SS_down_riboDbiased <- read_rds("../../output_data/aRMS_SS/aRMS_SS_down_riboDbiased.rds") # histogram
WT_NB_down_riboDbiased <- read_rds("../../output_data/NB_WT/WT_NB_down_riboDbiased.rds") # histogram
ALL_AML_down_riboDbiased <- read_rds("../../output_data/ALL_AML/ALL_AML_down_riboDbiased.rds") # histogram
```

### Fig S8

``` r
# for fig S8
# aRMS v SS
aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nSS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:28] $x
    INFO [2026-07-06 18:43:28] aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $category.names
    INFO [2026-07-06 18:43:28] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $filename
    INFO [2026-07-06 18:43:28] NULL
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $output
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $print.mode
    INFO [2026-07-06 18:43:28] c("raw", "percent")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fill
    INFO [2026-07-06 18:43:28] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.dist
    INFO [2026-07-06 18:43:28] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.default.pos
    INFO [2026-07-06 18:43:28] [1] "outer"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.pos
    INFO [2026-07-06 18:43:28] c(-14, 14, 175)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main
    INFO [2026-07-06 18:43:28] [1] "Upregulated genes in\nSS relative to aRMS"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.cex
    INFO [2026-07-06 18:43:28] [1] 1.1
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.pos
    INFO [2026-07-06 18:43:28] c(0.5, 1.2)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $disable.logging
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $sigdigs
    INFO [2026-07-06 18:43:28] [1] 3
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 

``` r
# NB v WT
WT_NB_up_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = WT_NB_up_polyAunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nWT relative to NB",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:28] $x
    INFO [2026-07-06 18:43:28] WT_NB_up_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $category.names
    INFO [2026-07-06 18:43:28] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $filename
    INFO [2026-07-06 18:43:28] NULL
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $output
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $print.mode
    INFO [2026-07-06 18:43:28] c("raw", "percent")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fill
    INFO [2026-07-06 18:43:28] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.dist
    INFO [2026-07-06 18:43:28] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.default.pos
    INFO [2026-07-06 18:43:28] [1] "outer"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.pos
    INFO [2026-07-06 18:43:28] c(-14, 14, 175)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main
    INFO [2026-07-06 18:43:28] [1] "Upregulated genes in\nWT relative to NB"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.cex
    INFO [2026-07-06 18:43:28] [1] 1.1
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.pos
    INFO [2026-07-06 18:43:28] c(0.5, 1.2)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $disable.logging
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $sigdigs
    INFO [2026-07-06 18:43:28] [1] 3
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 

``` r
# NB v WT
ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nAML relative to ALL",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:28] $x
    INFO [2026-07-06 18:43:28] ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $category.names
    INFO [2026-07-06 18:43:28] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $filename
    INFO [2026-07-06 18:43:28] NULL
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $output
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $print.mode
    INFO [2026-07-06 18:43:28] c("raw", "percent")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fill
    INFO [2026-07-06 18:43:28] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.cex
    INFO [2026-07-06 18:43:28] [1] 0.9
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.dist
    INFO [2026-07-06 18:43:28] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.default.pos
    INFO [2026-07-06 18:43:28] [1] "outer"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.pos
    INFO [2026-07-06 18:43:28] c(-14, 14, 175)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main
    INFO [2026-07-06 18:43:28] [1] "Upregulated genes in\nAML relative to ALL"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.cex
    INFO [2026-07-06 18:43:28] [1] 1.1
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.pos
    INFO [2026-07-06 18:43:28] c(0.5, 1.2)
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $disable.logging
    INFO [2026-07-06 18:43:28] [1] TRUE
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $sigdigs
    INFO [2026-07-06 18:43:28] [1] 3
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $cat.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $main.fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 
    INFO [2026-07-06 18:43:28] $fontfamily
    INFO [2026-07-06 18:43:28] [1] "sans"
    INFO [2026-07-06 18:43:28] 

``` r
Fig_S8 <- plot_grid(
  aRMS_SS_up_polyAunbiased_polyAbiased_riboDbiased_VD,
  WT_NB_up_polyAunbiased_polyAbiased_riboDbiased_VD,
  ALL_AML_up_polyAunbiased_polyAbiased_riboDbiased_VD,
  ncol = 3,
  label_size = 15)

ggsave("../../Figures/Fig_S8.png", Fig_S8, width = 10, height = 3, dpi = 300)
```

### Fig S9

``` r
# for fig S9
# aRMS v SS
aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nSS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:31] $x
    INFO [2026-07-06 18:43:31] aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $category.names
    INFO [2026-07-06 18:43:31] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $filename
    INFO [2026-07-06 18:43:31] NULL
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $output
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $print.mode
    INFO [2026-07-06 18:43:31] c("raw", "percent")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fill
    INFO [2026-07-06 18:43:31] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.dist
    INFO [2026-07-06 18:43:31] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.default.pos
    INFO [2026-07-06 18:43:31] [1] "outer"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.pos
    INFO [2026-07-06 18:43:31] c(-14, 14, 175)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main
    INFO [2026-07-06 18:43:31] [1] "Upregulated genes in\nSS relative to aRMS"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.cex
    INFO [2026-07-06 18:43:31] [1] 1.1
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.pos
    INFO [2026-07-06 18:43:31] c(0.5, 1.2)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $disable.logging
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $sigdigs
    INFO [2026-07-06 18:43:31] [1] 3
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 

``` r
# NB v WT
WT_NB_up_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = WT_NB_up_riboDunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nWT relative to NB",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:31] $x
    INFO [2026-07-06 18:43:31] WT_NB_up_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $category.names
    INFO [2026-07-06 18:43:31] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $filename
    INFO [2026-07-06 18:43:31] NULL
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $output
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $print.mode
    INFO [2026-07-06 18:43:31] c("raw", "percent")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fill
    INFO [2026-07-06 18:43:31] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.dist
    INFO [2026-07-06 18:43:31] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.default.pos
    INFO [2026-07-06 18:43:31] [1] "outer"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.pos
    INFO [2026-07-06 18:43:31] c(-14, 14, 175)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main
    INFO [2026-07-06 18:43:31] [1] "Upregulated genes in\nWT relative to NB"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.cex
    INFO [2026-07-06 18:43:31] [1] 1.1
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.pos
    INFO [2026-07-06 18:43:31] c(0.5, 1.2)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $disable.logging
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $sigdigs
    INFO [2026-07-06 18:43:31] [1] 3
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 

``` r
# NB v WT
ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Upregulated genes in\nAML relative to ALL",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:31] $x
    INFO [2026-07-06 18:43:31] ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $category.names
    INFO [2026-07-06 18:43:31] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $filename
    INFO [2026-07-06 18:43:31] NULL
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $output
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $print.mode
    INFO [2026-07-06 18:43:31] c("raw", "percent")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fill
    INFO [2026-07-06 18:43:31] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.cex
    INFO [2026-07-06 18:43:31] [1] 0.9
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.dist
    INFO [2026-07-06 18:43:31] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.default.pos
    INFO [2026-07-06 18:43:31] [1] "outer"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.pos
    INFO [2026-07-06 18:43:31] c(-14, 14, 175)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main
    INFO [2026-07-06 18:43:31] [1] "Upregulated genes in\nAML relative to ALL"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.cex
    INFO [2026-07-06 18:43:31] [1] 1.1
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.pos
    INFO [2026-07-06 18:43:31] c(0.5, 1.2)
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $disable.logging
    INFO [2026-07-06 18:43:31] [1] TRUE
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $sigdigs
    INFO [2026-07-06 18:43:31] [1] 3
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $cat.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $main.fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 
    INFO [2026-07-06 18:43:31] $fontfamily
    INFO [2026-07-06 18:43:31] [1] "sans"
    INFO [2026-07-06 18:43:31] 

``` r
Fig_S9 <- plot_grid(
  aRMS_SS_up_riboDunbiased_polyAbiased_riboDbiased_VD,
  WT_NB_up_riboDunbiased_polyAbiased_riboDbiased_VD,
  ALL_AML_up_riboDunbiased_polyAbiased_riboDbiased_VD,
  ncol = 3,
  label_size = 15)

ggsave("../../Figures/Fig_S9.png", Fig_S9, width = 10, height = 3, dpi = 300)
```

### Fig S10

``` r
# aRMS v SS
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
# NB v WT
WT_NB_up_polyAbiased_median <- WT_NB_up_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
WT_NB_up_polyAbiased_median
```

| median_value |
|-------------:|
|     1.430398 |

``` r
# AML v ALL
ALL_AML_up_polyAbiased_median <- ALL_AML_up_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
ALL_AML_up_polyAbiased_median
```

| median_value |
|-------------:|
|     1.342516 |

``` r
# aRMS v SS
aRMS_SS_up_polyAbiased_hist <- ggplot(aRMS_SS_up_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\naRMS v SS",
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

# NB v WT
WT_NB_up_polyAbiased_hist <- ggplot(WT_NB_up_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\nNB v WT",
#    subtitle = "polyA biased, Median LFC = 1.32",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = WT_NB_up_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=7, y=1080, label= "Median = 1.43")

# AML v ALL
ALL_AML_up_polyAbiased_hist <- ggplot(ALL_AML_up_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\nALL v AML",
#    subtitle = "polyA biased, Median LFC = 1.32",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = ALL_AML_up_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=3.7, y=580, label= "Median = 1.34")

Fig_S10 <- wrap_plots(
  aRMS_SS_up_polyAbiased_hist,
  WT_NB_up_polyAbiased_hist,
  ALL_AML_up_polyAbiased_hist,
  ncol = 3) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 4, r = 4),
    plot.title = element_text(hjust = 0.35, face = "bold", size = 14),
    legend.title = element_text(size = 12)
  )
Fig_S10
```

![](Supplemental_Figs_files/figure-commonmark/Fig_S10-1.png)

``` r
ggsave("../../Figures/Fig_S10.png", Fig_S10, width = 10, height = 3.4, dpi = 300)
```

### Fig S11

``` r
# aRMS v SS
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
# NB v WT
WT_NB_up_riboDbiased_median <- WT_NB_up_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
WT_NB_up_riboDbiased_median
```

| median_value |
|-------------:|
|     2.596998 |

``` r
# AML v ALL
ALL_AML_up_riboDbiased_median <- ALL_AML_up_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
ALL_AML_up_riboDbiased_median
```

| median_value |
|-------------:|
|     1.752502 |

``` r
# aRMS v SS
aRMS_SS_up_riboDbiased_hist <- ggplot(aRMS_SS_up_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\naRMS v SS",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = aRMS_SS_up_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=8, y=750, label= "Median = 2.76")


# NB v WT
WT_NB_up_riboDbiased_hist <- ggplot(WT_NB_up_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\nNB v WT",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = WT_NB_up_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=6, y=700, label= "Median = 2.60")

# AML v ALL
ALL_AML_up_riboDbiased_hist <- ggplot(ALL_AML_up_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Upregulated in\nAML v ALL",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = ALL_AML_up_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=5, y=500, label= "Median = 1.75")

Fig_S11 <- wrap_plots(
  aRMS_SS_up_riboDbiased_hist,
  WT_NB_up_riboDbiased_hist,
  ALL_AML_up_riboDbiased_hist,
  ncol = 3) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 4, r = 4),
    plot.title = element_text(hjust = 0.35, face = "bold", size = 14),
    legend.title = element_text(size = 12)
  )
Fig_S11
```

![](Supplemental_Figs_files/figure-commonmark/Fig_S11-1.png)

``` r
ggsave("../../Figures/Fig_S11.png", Fig_S11, width = 10, height = 3.4, dpi = 300)
```

### Fig S12

``` r
# for fig S12 - downregulated genes
# aRMS v SS
aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nSS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:38] $x
    INFO [2026-07-06 18:43:38] aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $category.names
    INFO [2026-07-06 18:43:38] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $filename
    INFO [2026-07-06 18:43:38] NULL
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $output
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $print.mode
    INFO [2026-07-06 18:43:38] c("raw", "percent")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fill
    INFO [2026-07-06 18:43:38] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.dist
    INFO [2026-07-06 18:43:38] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.default.pos
    INFO [2026-07-06 18:43:38] [1] "outer"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.pos
    INFO [2026-07-06 18:43:38] c(-14, 14, 175)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main
    INFO [2026-07-06 18:43:38] [1] "Downregulated genes in\nSS relative to aRMS"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.cex
    INFO [2026-07-06 18:43:38] [1] 1.1
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.pos
    INFO [2026-07-06 18:43:38] c(0.5, 1.2)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $disable.logging
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $sigdigs
    INFO [2026-07-06 18:43:38] [1] 3
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 

``` r
# NB v WT
WT_NB_down_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = WT_NB_down_polyAunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nWT relative to NB",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:38] $x
    INFO [2026-07-06 18:43:38] WT_NB_down_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $category.names
    INFO [2026-07-06 18:43:38] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $filename
    INFO [2026-07-06 18:43:38] NULL
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $output
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $print.mode
    INFO [2026-07-06 18:43:38] c("raw", "percent")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fill
    INFO [2026-07-06 18:43:38] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.dist
    INFO [2026-07-06 18:43:38] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.default.pos
    INFO [2026-07-06 18:43:38] [1] "outer"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.pos
    INFO [2026-07-06 18:43:38] c(-14, 14, 175)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main
    INFO [2026-07-06 18:43:38] [1] "Downregulated genes in\nWT relative to NB"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.cex
    INFO [2026-07-06 18:43:38] [1] 1.1
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.pos
    INFO [2026-07-06 18:43:38] c(0.5, 1.2)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $disable.logging
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $sigdigs
    INFO [2026-07-06 18:43:38] [1] 3
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 

``` r
# NB v WT
ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nAML relative to ALL",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:38] $x
    INFO [2026-07-06 18:43:38] ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $category.names
    INFO [2026-07-06 18:43:38] c("PolyA unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $filename
    INFO [2026-07-06 18:43:38] NULL
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $output
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $print.mode
    INFO [2026-07-06 18:43:38] c("raw", "percent")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fill
    INFO [2026-07-06 18:43:38] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.cex
    INFO [2026-07-06 18:43:38] [1] 0.9
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.dist
    INFO [2026-07-06 18:43:38] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.default.pos
    INFO [2026-07-06 18:43:38] [1] "outer"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.pos
    INFO [2026-07-06 18:43:38] c(-14, 14, 175)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main
    INFO [2026-07-06 18:43:38] [1] "Downregulated genes in\nAML relative to ALL"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.cex
    INFO [2026-07-06 18:43:38] [1] 1.1
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.pos
    INFO [2026-07-06 18:43:38] c(0.5, 1.2)
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $disable.logging
    INFO [2026-07-06 18:43:38] [1] TRUE
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $sigdigs
    INFO [2026-07-06 18:43:38] [1] 3
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $cat.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $main.fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 
    INFO [2026-07-06 18:43:38] $fontfamily
    INFO [2026-07-06 18:43:38] [1] "sans"
    INFO [2026-07-06 18:43:38] 

``` r
Fig_S12 <- plot_grid(
  aRMS_SS_down_polyAunbiased_polyAbiased_riboDbiased_VD,
  WT_NB_down_polyAunbiased_polyAbiased_riboDbiased_VD,
  ALL_AML_down_polyAunbiased_polyAbiased_riboDbiased_VD,
  ncol = 3,
  label_size = 15)

ggsave("../../Figures/Fig_S12.png", Fig_S12, width = 10, height = 3, dpi = 300)
```

### Fig S13

``` r
# for fig S13
# aRMS v SS
aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nSS relative to aRMS",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:40] $x
    INFO [2026-07-06 18:43:40] aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $category.names
    INFO [2026-07-06 18:43:40] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $filename
    INFO [2026-07-06 18:43:40] NULL
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $output
    INFO [2026-07-06 18:43:40] [1] TRUE
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $print.mode
    INFO [2026-07-06 18:43:40] c("raw", "percent")
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $fill
    INFO [2026-07-06 18:43:40] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cex
    INFO [2026-07-06 18:43:40] [1] 0.9
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cat.cex
    INFO [2026-07-06 18:43:40] [1] 0.9
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cat.dist
    INFO [2026-07-06 18:43:40] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cat.default.pos
    INFO [2026-07-06 18:43:40] [1] "outer"
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cat.pos
    INFO [2026-07-06 18:43:40] c(-14, 14, 175)
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $main
    INFO [2026-07-06 18:43:40] [1] "Downregulated genes in\nSS relative to aRMS"
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $main.cex
    INFO [2026-07-06 18:43:40] [1] 1.1
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $main.pos
    INFO [2026-07-06 18:43:40] c(0.5, 1.2)
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $disable.logging
    INFO [2026-07-06 18:43:40] [1] TRUE
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $sigdigs
    INFO [2026-07-06 18:43:40] [1] 3
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $cat.fontfamily
    INFO [2026-07-06 18:43:40] [1] "sans"
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $main.fontfamily
    INFO [2026-07-06 18:43:40] [1] "sans"
    INFO [2026-07-06 18:43:40] 
    INFO [2026-07-06 18:43:40] $fontfamily
    INFO [2026-07-06 18:43:40] [1] "sans"
    INFO [2026-07-06 18:43:40] 

``` r
# NB v WT
WT_NB_down_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = WT_NB_down_riboDunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nWT relative to NB",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:41] $x
    INFO [2026-07-06 18:43:41] WT_NB_down_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $category.names
    INFO [2026-07-06 18:43:41] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $filename
    INFO [2026-07-06 18:43:41] NULL
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $output
    INFO [2026-07-06 18:43:41] [1] TRUE
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $print.mode
    INFO [2026-07-06 18:43:41] c("raw", "percent")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $fill
    INFO [2026-07-06 18:43:41] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cex
    INFO [2026-07-06 18:43:41] [1] 0.9
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.cex
    INFO [2026-07-06 18:43:41] [1] 0.9
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.dist
    INFO [2026-07-06 18:43:41] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.default.pos
    INFO [2026-07-06 18:43:41] [1] "outer"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.pos
    INFO [2026-07-06 18:43:41] c(-14, 14, 175)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main
    INFO [2026-07-06 18:43:41] [1] "Downregulated genes in\nWT relative to NB"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.cex
    INFO [2026-07-06 18:43:41] [1] 1.1
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.pos
    INFO [2026-07-06 18:43:41] c(0.5, 1.2)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $disable.logging
    INFO [2026-07-06 18:43:41] [1] TRUE
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $sigdigs
    INFO [2026-07-06 18:43:41] [1] 3
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 

``` r
# NB v WT
ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased_VD <- venn.diagram(
  x = ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased,
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
  cex = 0.9, # Font size for counts
  cat.cex = 0.9, # Font size for category names
  cat.dist = c(0.04, 0.04, 0.04),
  # height = 1000,
  # width = 1000,
  cat.default.pos = "outer",
  cat.pos = c(-14, 14, 175), 
  main = "Downregulated genes in\nAML relative to ALL",
  main.cex = 1.1, # Font size for main title
  main.pos = c(0.5, 1.2), 
  disable.logging = TRUE,
  sigdigs = 3, 
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  fontfamily = "sans"
)
```

    INFO [2026-07-06 18:43:41] $x
    INFO [2026-07-06 18:43:41] ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $category.names
    INFO [2026-07-06 18:43:41] c("RiboD unbiased", "PolyA biased", "RiboD biased")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $filename
    INFO [2026-07-06 18:43:41] NULL
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $output
    INFO [2026-07-06 18:43:41] [1] TRUE
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $print.mode
    INFO [2026-07-06 18:43:41] c("raw", "percent")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $fill
    INFO [2026-07-06 18:43:41] c("#BB5566", "#0072B2", "#E69F00")
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cex
    INFO [2026-07-06 18:43:41] [1] 0.9
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.cex
    INFO [2026-07-06 18:43:41] [1] 0.9
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.dist
    INFO [2026-07-06 18:43:41] c(0.04, 0.04, 0.04)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.default.pos
    INFO [2026-07-06 18:43:41] [1] "outer"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.pos
    INFO [2026-07-06 18:43:41] c(-14, 14, 175)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main
    INFO [2026-07-06 18:43:41] [1] "Downregulated genes in\nAML relative to ALL"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.cex
    INFO [2026-07-06 18:43:41] [1] 1.1
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.pos
    INFO [2026-07-06 18:43:41] c(0.5, 1.2)
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $disable.logging
    INFO [2026-07-06 18:43:41] [1] TRUE
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $sigdigs
    INFO [2026-07-06 18:43:41] [1] 3
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $cat.fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $main.fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 
    INFO [2026-07-06 18:43:41] $fontfamily
    INFO [2026-07-06 18:43:41] [1] "sans"
    INFO [2026-07-06 18:43:41] 

``` r
Fig_S13 <- plot_grid(
  aRMS_SS_down_riboDunbiased_polyAbiased_riboDbiased_VD,
  WT_NB_down_riboDunbiased_polyAbiased_riboDbiased_VD,
  ALL_AML_down_riboDunbiased_polyAbiased_riboDbiased_VD,
  ncol = 3,
  label_size = 15)

ggsave("../../Figures/Fig_S13.png", Fig_S13, width = 10, height = 3, dpi = 300)
```

### Fig S14

``` r
# aRMS v SS
aRMS_SS_down_polyAbiased_median <- aRMS_SS_down_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
aRMS_SS_down_polyAbiased_median
```

| median_value |
|-------------:|
|    -2.895363 |

``` r
# NB v WT
WT_NB_down_polyAbiased_median <- WT_NB_down_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
WT_NB_down_polyAbiased_median
```

| median_value |
|-------------:|
|    -3.066331 |

``` r
# AML v ALL
ALL_AML_down_polyAbiased_median <- ALL_AML_down_polyAbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
ALL_AML_down_polyAbiased_median
```

| median_value |
|-------------:|
|    -1.627678 |

``` r
# aRMS v SS
aRMS_SS_down_polyAbiased_hist <- ggplot(aRMS_SS_down_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\naRMS v SS",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = aRMS_SS_down_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-7, y=580, label= "Median = -2.9")

# NB v WT
WT_NB_down_polyAbiased_hist <- ggplot(WT_NB_down_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\nNB v WT",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = WT_NB_down_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-10, y=750, label= "Median = -3.1")

# AML v ALL
ALL_AML_down_polyAbiased_hist <- ggplot(ALL_AML_down_polyAbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#0072B2",
    color = "#0072B2",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\nALL v AML",
    y = "Number of genes unique to PolyA biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = ALL_AML_down_polyAbiased_median$median_value), color = "#004166", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-6.7, y=750, label= "Median = -1.62")

Fig_S14 <- wrap_plots(
  aRMS_SS_down_polyAbiased_hist,
  WT_NB_down_polyAbiased_hist,
  ALL_AML_down_polyAbiased_hist,
  ncol = 3) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 4, r = 4),
    plot.title = element_text(hjust = 0.35, face = "bold", size = 14),
    legend.title = element_text(size = 12)
  )
Fig_S14
```

![](Supplemental_Figs_files/figure-commonmark/Fig_S14-1.png)

``` r
ggsave("../../Figures/Fig_S14.png", Fig_S14, width = 10, height = 3.4, dpi = 300)
```

### Fig S15

``` r
# aRMS v SS
aRMS_SS_down_riboDbiased_median <- aRMS_SS_down_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
aRMS_SS_down_riboDbiased_median
```

| median_value |
|-------------:|
|    -1.355701 |

``` r
# NB v WT
WT_NB_down_riboDbiased_median <- WT_NB_down_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
WT_NB_down_riboDbiased_median
```

| median_value |
|-------------:|
|    -1.334511 |

``` r
# AML v ALL
ALL_AML_down_riboDbiased_median <- ALL_AML_down_riboDbiased %>%
  summarise(
    median_value = median(log2FoldChange.x),
  )
ALL_AML_down_riboDbiased_median
```

| median_value |
|-------------:|
|    -1.296047 |

``` r
# aRMS v SS
aRMS_SS_down_riboDbiased_hist <- ggplot(aRMS_SS_down_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\naRMS v SS",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = aRMS_SS_down_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-8, y=700, label= "Median = -1.36")


# NB v WT
WT_NB_down_riboDbiased_hist <- ggplot(WT_NB_down_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\nNB v WT",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = WT_NB_down_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-5, y=400, label= "Median = -1.33")

# AML v ALL
ALL_AML_down_riboDbiased_hist <- ggplot(ALL_AML_down_riboDbiased, aes(x = log2FoldChange.x)) +
  geom_histogram(
    bins = 30, # 30 is default
    fill = "#E69F00",
    color = "#E69F00",
    alpha = 0.7
    ) +
#  scale_y_continuous(breaks = c(200, 400, 600, 800), limits = c(0,1000)) +
#  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16)) +
  labs(
    title = "Downregulated in\nAML v ALL",
    y = "Number of genes unique to RiboD biased",
    x = "log2 Fold Change"
) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
        ) +
  geom_vline(aes(xintercept = ALL_AML_down_riboDbiased_median$median_value), color = "#9a6a00", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x=-5, y=400, label= "Median = -1.30")

Fig_S15 <- wrap_plots(
  aRMS_SS_down_riboDbiased_hist,
  WT_NB_down_riboDbiased_hist,
  ALL_AML_down_riboDbiased_hist,
  ncol = 3) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") &
  theme(
    plot.margin = margin(t = 2, b = 2, l = 4, r = 4),
    plot.title = element_text(hjust = 0.35, face = "bold", size = 14),
    legend.title = element_text(size = 12)
  )
Fig_S15
```

![](Supplemental_Figs_files/figure-commonmark/Fig_S15-1.png)

``` r
ggsave("../../Figures/Fig_S15.png", Fig_S15, width = 10, height = 3.4, dpi = 300)
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
     date     2026-07-06
     pandoc   3.8.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
     quarto   1.9.36 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     ! package        * version date (UTC) lib source
     P cli              3.6.5   2025-04-23 [?] CRAN (R 4.5.0)
     P cowplot        * 1.2.0   2025-07-07 [?] RSPM
     P digest           0.6.37  2024-08-19 [?] CRAN (R 4.5.0)
     P dplyr          * 1.1.4   2023-11-17 [?] CRAN (R 4.5.0)
     P evaluate         1.0.5   2025-08-27 [?] RSPM
     P farver           2.1.2   2024-05-13 [?] RSPM
     P fastmap          1.2.0   2024-05-15 [?] RSPM
     P forcats        * 1.0.1   2025-09-25 [?] RSPM
     P formatR          1.14    2023-01-17 [?] CRAN (R 4.5.0)
     P futile.logger  * 1.4.3   2016-07-10 [?] CRAN (R 4.5.0)
     P futile.options   1.0.1   2018-04-20 [?] CRAN (R 4.5.0)
     P generics         0.1.4   2025-05-09 [?] RSPM
     P ggplot2        * 4.0.0   2025-09-11 [?] CRAN (R 4.5.0)
     P glue             1.8.0   2024-09-30 [?] CRAN (R 4.5.0)
     P gtable           0.3.6   2024-10-25 [?] RSPM
     P hms              1.1.4   2025-10-17 [?] RSPM
     P htmltools        0.5.8.1 2024-04-04 [?] CRAN (R 4.5.0)
     P jsonlite         2.0.0   2025-03-27 [?] RSPM
     P knitr            1.50    2025-03-16 [?] CRAN (R 4.5.0)
     P labeling         0.4.3   2023-08-29 [?] RSPM
     P lambda.r         1.2.4   2019-09-18 [?] CRAN (R 4.5.0)
     P lifecycle        1.0.4   2023-11-07 [?] CRAN (R 4.5.0)
     P lubridate      * 1.9.4   2024-12-08 [?] CRAN (R 4.5.0)
     P magrittr         2.0.4   2025-09-12 [?] CRAN (R 4.5.0)
     P patchwork      * 1.3.2   2025-08-25 [?] RSPM
     P pillar           1.11.1  2025-09-17 [?] RSPM
     P pkgconfig        2.0.3   2019-09-22 [?] RSPM
     P purrr          * 1.1.0   2025-07-10 [?] CRAN (R 4.5.0)
     P R6               2.6.1   2025-02-15 [?] RSPM
     P ragg             1.5.0   2025-09-02 [?] CRAN (R 4.5.0)
     P RColorBrewer     1.1-3   2022-04-03 [?] RSPM
     P readr          * 2.1.5   2024-01-10 [?] CRAN (R 4.5.0)
       renv             1.1.5   2025-07-24 [1] CRAN (R 4.5.0)
     P rlang            1.2.0   2026-04-06 [?] RSPM
     P rmarkdown        2.30    2025-09-28 [?] CRAN (R 4.5.0)
     P rstudioapi       0.17.1  2024-10-22 [?] CRAN (R 4.5.0)
     P S7               0.2.0   2024-11-07 [?] CRAN (R 4.5.0)
     P scales           1.4.0   2025-04-24 [?] RSPM
     P sessioninfo      1.2.3   2025-02-05 [?] CRAN (R 4.5.0)
     P stringi          1.8.7   2025-03-27 [?] RSPM
     P stringr        * 1.5.2   2025-09-08 [?] CRAN (R 4.5.0)
     P systemfonts      1.3.1   2025-10-01 [?] CRAN (R 4.5.0)
     P textshaping      1.0.4   2025-10-10 [?] CRAN (R 4.5.0)
     P tibble         * 3.3.0   2025-06-08 [?] CRAN (R 4.5.0)
     P tidyr          * 1.3.1   2024-01-24 [?] CRAN (R 4.5.0)
     P tidyselect       1.2.1   2024-03-11 [?] RSPM
     P tidyverse      * 2.0.0   2023-02-22 [?] RSPM
     P timechange       0.3.0   2024-01-18 [?] CRAN (R 4.5.0)
     P tzdb             0.5.0   2025-03-15 [?] RSPM
     P vctrs            0.6.5   2023-12-01 [?] CRAN (R 4.5.0)
     P VennDiagram    * 1.8.2   2026-01-11 [?] RSPM
     P withr            3.0.2   2024-10-28 [?] CRAN (R 4.5.0)
     P xfun             0.55    2025-12-16 [?] CRAN (R 4.5.2)
     P yaml             2.3.10  2024-07-26 [?] CRAN (R 4.5.0)

     [1] /Users/maryke/Documents/Treehouse/Lab_Notebooks/transcript_enrichment_bias_assessment/Fig_2/Supplemental_Figs/renv/library/macos/R-4.5/aarch64-apple-darwin20
     [2] /Users/maryke/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74

     * ── Packages attached to the search path.
     P ── Loaded and on-disk path mismatch.

    ──────────────────────────────────────────────────────────────────────────────
