# assemble-expression-for-THR13 0970-polyA-vs-riboD-comparison-v2
Holly Beale
2026-09-16

- [Import and groom gene expression
  levels](#import-and-groom-gene-expression-levels)

``` r
gene_names <- read_tsv("EnsGeneID_Hugo_Observed_Conversions.txt.gz", show_col_types = FALSE)

treehouse_druggable_genes <- read_tsv("treehouseDruggableGenes_2020-03_25.txt")
```

    Rows: 115 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): gene, group

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
meta <- read_tsv("UCSC in-house RNA-seq validation samples - assign thids.tsv")
```

    Rows: 12 Columns: 5
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (3): th_sampleid, File name, Library Slxn Method
    dbl (2): donor_key, set row background color (formatting purposes only!)

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
meta %>%
  select(th_sampleid, `Library Slxn Method`, `File name`) %>%
  filter(str_detect(`File name`, "_R1_")) 
```

| th_sampleid    | Library Slxn Method | File name                          |
|:---------------|:--------------------|:-----------------------------------|
| THR13_0970_S12 | polyA               | DIPG4_polyA_001_S1_R1_001.fastq.gz |
| THR13_0970_S13 | polyA               | DIPG4_polyA_002_S2_R1_001.fastq.gz |
| THR13_0970_S14 | polyA               | DIPG4_polyA_003_S3_R1_001.fastq.gz |
| THR13_0970_S15 | riboD               | DIPG4_riboD_001_S4_R1_001.fastq.gz |
| THR13_0970_S16 | riboD               | DIPG4_riboD_002_S5_R1_001.fastq.gz |
| THR13_0970_S17 | riboD               | DIPG4_riboD_003_S6_R1_001.fastq.gz |

# Import and groom gene expression levels

``` r
rsem_pattern <- "_rsem_genes.results"
rsem_files <- list.files(pattern = rsem_pattern)
names(rsem_files) <- str_remove(rsem_files, rsem_pattern)

rsem_raw_data <- rsem_files %>%
  purrr::map_dfr(read_tsv, show_col_types = FALSE, .id="sample") %>%
  select(-`transcript_id(s)`)

rsem_ens_log2TPM1 <- rsem_raw_data %>%
  left_join(gene_names, by=c("gene_id"= "EnsGeneID")) %>%
  mutate(log2TPM1 = log2(TPM+1)) %>%
  left_join(meta %>% 
              select(th_sampleid, lib_prep = `Library Slxn Method`) %>%
              distinct(),
            by=c("sample" = "th_sampleid")) %>%
  rename(full_sample = sample) %>%
  mutate(sample = str_remove(full_sample, "THR13_0970_"),
         is_treehouse_druggable_gene = HugoID %in% treehouse_druggable_genes$gene) %>% 
  select(full_sample, ensembl_gene_ID = gene_id, HugoID,  TPM, 
         log2TPM1, lib_prep, sample, is_treehouse_druggable_gene)
```

``` r
write_tsv(rsem_ens_log2TPM1, "rsem_ensembl_log2TPM1_THR13_0970_S12-S17.tsv.gz")
```

``` r
nrow(rsem_ens_log2TPM1)/6
```

    [1] 60498
