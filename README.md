# transcript_enrichment_bias_assessment
This is a comparative assessment of PolyA and RiboD transcript enrichment methods on gene expression profiles in pediatric tumor samples from the Treehouse v25.01 Tumor Compendia. This 
This code generates the figures and data found in [ link to paper ], with the exception of the UMAPs in figure 1, which is found at https://github.com/UCSC-Treehouse/lib-prep-visualization . The abstract for the manuscript is as follows:

## Background
RNA-sequencing has proven to be an invaluable tool for biomedical research, but transcript enrichment methods have yet to be standardized in the clinical setting. Two widely used methods include polyadenylation-selection (PolyA) and ribosomal RNA-depletion (RiboD). These methods collect different populations of RNA, thus limiting comparable samples to those prepared with the same method. Therefore, a clear understanding of the transcript enrichment method used in the library preparation is essential to comparative gene expression analyses. In this study, we evaluated the effects of transcript enrichment methods on measurable gene expression and quantified the biases introduced when comparing samples prepared with different transcript enrichment methods. Using the Treehouse Childhood Cancer Initiative’s publicly available bulk RNA-sequencing compendia, we analyzed the effects of different transcript enrichment methods on sample similarity, gene expression profiles, and differential expression analysis.
## Findings
We found that samples prepared from the same transcript enrichment method are similar within a given disease, but samples form distinct groupings based on transcript enrichment method when comparing across methods. PolyA and RiboD samples yield different amounts of measurable genes within a given disease, with variability in gene expression between PolyA and RiboD samples for a given gene.  These different transcriptomic profiles result in conflicting numbers of differentially expressed genes when comparing diseases across transcript enrichment methods.
## Conclusions
Our study indicates that samples prepared with different transcript enrichment methods yield distinct transcriptomic profiles. Consequently, caution must be taken when performing comparative gene expression studies across disease types, and the usage of transcript enrichment methods must be well-documented and considered prior to comparative analysis.

## Directory structure
```
transcript_enrichment_bias_assessment/
├── .gitignore
├── Fig_1/                  # contains scripts used for analysis in Fig1
├── Fig_2/                  # contains scripts used for DESeq2 analysis in Fig2
├── LICENSE
├── README.md
├── correlation_analysis/   # contains initial correlation analysis script for sample selection
└── output_data/
    ├── ALL_AML/            # contains DESeq2 output data for ALL vs AML
    ├── Fig1C/
    ├── Fig1D/
    ├── Fig1E_F/
    ├── Fig1G/
    ├── Fig1H/
    ├── NB_WT/              # contains DESeq2 output data for NB vs WT
    ├── aRMS_SS/            # contains DESeq2 output data for aRMS vs SS
    └── all_disease_comparisons/  # contains uniquely up/downregulated gene lists from DESeq2 analysis across disease comparisons
```
