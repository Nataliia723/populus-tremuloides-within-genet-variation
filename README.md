# Within-genet genotypic variation in trembling aspen (*Populus tremuloides*)

This repository contains the R code used to investigate within-genet genotypic variation in trembling aspen (*Populus tremuloides*), as reported in the manuscript and its Supporting Information.

## Analyses

The R script performs the following analyses:

- Exploratory principal component analyses (PCA) within genets
- Identification of genotype discordances relative to within-genet consensus genotypes
- Quality filtering based on sequencing depth (DP), genotype quality (GQ), allele depth (AD), and allele balance
- Identification of recurrent high-quality genotype discordances
- Detection of local multilocus genotypic signatures
- Generation of summary tables and supplementary results

The analyses examine genotypic variation within identified genets. Detected discordances are not interpreted as confirmed somatic mutations.

## R script

The main analysis script is:

`within_genet_genotypic_variation.R`

## Dataset

The analyses were conducted using:

- 421 trembling aspen samples
- 48,976 unlinked SNPs
- 33 identified genets

Within-genet analyses were performed on genets represented by at least three samples.

## Software and R packages

The analyses were performed using R and require the following packages:

- dplyr
- tidyr
- ggplot2
- purrr
- readr

The packages can be installed in R using:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "readr"
))
```

## Input files

The analyses require the following input data:

**Genomic dataset**

`second_filters_m4_p80_x0_s3.canonical.unlinked_0.5_100k.vcf`

Available from Zenodo:

https://doi.org/10.5281/zenodo.23169870

**Genet assignments**

A table containing sample identifiers and their corresponding genet assignments.

The analysis script expects a data frame named `genet_map` containing the following columns:

- Sample
- Genet

**VCF data**

The script requires a data frame named `vcf` containing the genomic data and sample identifiers.

Genotype information includes:

- GT: Genotype
- DP: Sequencing depth
- AD: Allele depth
- GQ: Genotype quality

## Quality filtering

Genotype discordances were evaluated using the following criteria:

- Sequencing depth: DP ≥ 10
- Genotype quality: GQ ≥ 30
- Homozygous reference (0/0): allele balance ≤ 0.05
- Heterozygous (0/1): allele balance between 0.25 and 0.75
- Homozygous alternative (1/1): allele balance ≥ 0.95

Recurrent discordances were examined across samples within the same genet.

## Reproducibility

The R code corresponds to the within-genet genotypic variation analyses described in the manuscript and Supporting Information.

The script generates summary tables and output files documenting genotype discordances, quality filtering, recurrent loci, and local multilocus signatures.

The input objects `vcf` and `genet_map` must be loaded into the R environment before running the script.

## Code availability

The R code used for within-genet genotypic variation analyses is publicly available in this repository.
