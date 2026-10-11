# ================================================================
# WITHIN-GENET GENOTYPIC VARIATION IN POPULUS TREMULOIDES

# This script evaluates within-genet genotypic variation in
# Populus tremuloides using the unlinked dataset of 48,976 SNPs.
#
# Analyses include:
# - exploratory PCA within genets;
# - genotype discordance relative to the within-genet consensus;
# - quality filtering using DP, GQ, AD and allele balance;
# - recurrence of high-quality discordances among samples;
# - physical clustering and identification of recurrent local
#   multilocus genotypic signatures.
#
# High-quality discordances are interpreted as within-genet
# genotypic differences, not as confirmed somatic mutations.
#
# Dataset:
# 421 samples
# 48,976 unlinked SNPs
#
# VCF genotype information:
# GT = genotype (e.g., 0/0, 0/1, 1/1)
# DP = sequencing depth (total number of reads covering the locus)
# AD = allele depth (number of reads supporting the reference and alternative alleles)
# GQ = genotype quality (confidence in the assigned genotype)
# ================================================================


# ================================================================
# 0. PACKAGES
# ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)


# ================================================================
# 1. INPUT OBJECTS
# ================================================================

# Required objects:
#
# vcf
#   data.frame containing the unlinked VCF
#   48,976 SNPs × metadata + 421 sample columns
#
# genet_map
#   data.frame/tibble with:
#   Sample
#   Genet
#
# vcf_fields
#   list containing matrices:
#   $GT
#   $DP
#   $AD_REF
#   $AD_ALT
#   $GQ
#
# Sample names are expected to begin with "ptrem_"


# ================================================================
# 2. EXTRACT GT, DP, AD AND GQ FROM VCF
#    Run only if vcf_fields does not already exist
# ================================================================

extract_vcf_fields <- function(vcf) {
  
  sample_cols <- grep(
    "^ptrem_",
    names(vcf),
    value = TRUE
  )
  
  n_snp <- nrow(vcf)
  n_sample <- length(sample_cols)
  
  GT <- matrix(
    NA_character_,
    nrow = n_snp,
    ncol = n_sample,
    dimnames = list(NULL, sample_cols)
  )
  
  DP <- matrix(
    NA_real_,
    nrow = n_snp,
    ncol = n_sample,
    dimnames = list(NULL, sample_cols)
  )
  
  AD_REF <- matrix(
    NA_real_,
    nrow = n_snp,
    ncol = n_sample,
    dimnames = list(NULL, sample_cols)
  )
  
  AD_ALT <- matrix(
    NA_real_,
    nrow = n_snp,
    ncol = n_sample,
    dimnames = list(NULL, sample_cols)
  )
  
  GQ <- matrix(
    NA_real_,
    nrow = n_snp,
    ncol = n_sample,
    dimnames = list(NULL, sample_cols)
  )
  
  for (j in seq_along(sample_cols)) {
    
    z <- strsplit(
      as.character(vcf[[sample_cols[j]]]),
      ":",
      fixed = TRUE
    )
    
    GT[, j] <- vapply(
      z,
      function(x) {
        if (length(x) >= 1 && x[1] != "./.") x[1]
        else NA_character_
      },
      character(1)
    )
    
    DP[, j] <- suppressWarnings(
      as.numeric(
        vapply(
          z,
          function(x) {
            if (length(x) >= 2) x[2]
            else NA_character_
          },
          character(1)
        )
      )
    )
    
    AD_string <- vapply(
      z,
      function(x) {
        if (length(x) >= 3) x[3]
        else NA_character_
      },
      character(1)
    )
    
    AD_split <- strsplit(
      AD_string,
      ",",
      fixed = TRUE
    )
    
    AD_REF[, j] <- suppressWarnings(
      as.numeric(
        vapply(
          AD_split,
          function(x) {
            if (length(x) >= 1) x[1]
            else NA_character_
          },
          character(1)
        )
      )
    )
    
    AD_ALT[, j] <- suppressWarnings(
      as.numeric(
        vapply(
          AD_split,
          function(x) {
            if (length(x) >= 2) x[2]
            else NA_character_
          },
          character(1)
        )
      )
    )
    
    GQ[, j] <- suppressWarnings(
      as.numeric(
        vapply(
          z,
          function(x) {
            if (length(x) >= 4) x[4]
            else NA_character_
          },
          character(1)
        )
      )
    )
  }
  
  list(
    GT = GT,
    DP = DP,
    AD_REF = AD_REF,
    AD_ALT = AD_ALT,
    GQ = GQ
  )
}


if (!exists("vcf_fields")) {
  vcf_fields <- extract_vcf_fields(vcf)
}


# ================================================================
# 3. BASIC CHECKS
# ================================================================

dim(vcf_fields$GT)
dim(vcf_fields$DP)
dim(vcf_fields$AD_REF)
dim(vcf_fields$AD_ALT)
dim(vcf_fields$GQ)

head(genet_map)

stopifnot(
  identical(
    colnames(vcf_fields$GT),
    colnames(vcf_fields$DP)
  )
)

stopifnot(
  identical(
    colnames(vcf_fields$GT),
    colnames(vcf_fields$GQ)
  )
)


# ================================================================
# 4. CONVERT GT TO NUMERIC
#
# 0/0 = 0
# 0/1 = 1
# 1/0 = 1
# 1/1 = 2
# ================================================================

GT_num <- matrix(
  NA_real_,
  nrow = nrow(vcf_fields$GT),
  ncol = ncol(vcf_fields$GT),
  dimnames = dimnames(vcf_fields$GT)
)

GT_num[vcf_fields$GT == "0/0"] <- 0

GT_num[
  vcf_fields$GT %in% c("0/1", "1/0")
] <- 1

GT_num[vcf_fields$GT == "1/1"] <- 2


# ================================================================
# 5. CONSENSUS GENOTYPE FUNCTION
#
# Ties are returned as NA.
# ================================================================

consensus_gt <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  tab <- table(x)
  
  if (sum(tab == max(tab)) > 1) {
    return(NA_real_)
  }
  
  as.numeric(
    names(tab)[which.max(tab)]
  )
}


# ================================================================
# PART A
# EXPLORATORY WITHIN-GENET PCA
# ================================================================


# ================================================================
# 6. FUNCTION FOR WITHIN-GENET PCA
#
# Criteria:
# - genet must contain >= 3 samples
# - SNP call rate >= 80% within genet
# - retain variable SNPs
# - mean-impute remaining missing values
# - center SNPs
# - do NOT scale SNPs
#
# PCA is exploratory and does NOT identify somatic mutations.
# ================================================================

run_within_genet_pca <- function(
    genet_name,
    geno = GT_num,
    genet_map,
    min_call_rate = 0.80) {
  
  samples <- genet_map$Sample[
    genet_map$Genet == genet_name
  ]
  
  samples <- intersect(
    samples,
    colnames(geno)
  )
  
  if (length(samples) < 3) {
    return(NULL)
  }
  
  mat <- geno[
    ,
    samples,
    drop = FALSE
  ]
  
  # SNP call rate
  call_rate <- rowMeans(
    !is.na(mat)
  )
  
  mat <- mat[
    call_rate >= min_call_rate,
    ,
    drop = FALSE
  ]
  
  # Retain variable SNPs
  variable <- apply(
    mat,
    1,
    function(x) {
      length(unique(x[!is.na(x)])) > 1
    }
  )
  
  mat <- mat[
    variable,
    ,
    drop = FALSE
  ]
  
  if (nrow(mat) < 2) {
    return(NULL)
  }
  
  # Mean imputation by SNP
  for (i in seq_len(nrow(mat))) {
    
    missing <- is.na(mat[i, ])
    
    if (any(missing)) {
      
      mat[
        i,
        missing
      ] <- mean(
        mat[i, ],
        na.rm = TRUE
      )
    }
  }
  
  # prcomp expects samples in rows
  pca <- prcomp(
    t(mat),
    center = TRUE,
    scale. = FALSE
  )
  
  variance <- (
    pca$sdev^2 /
      sum(pca$sdev^2)
  ) * 100
  
  scores <- data.frame(
    Genet = genet_name,
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  
  list(
    Genet = genet_name,
    N_samples = length(samples),
    N_variable_SNPs = nrow(mat),
    PC1_percent = variance[1],
    PC2_percent = variance[2],
    scores = scores,
    pca = pca
  )
}


# ================================================================
# 7. RUN PCA FOR ALL GENETS WITH >= 3 SAMPLES
# ================================================================

genets_pca <- genet_map %>%
  count(
    Genet,
    name = "N_samples"
  ) %>%
  filter(
    N_samples >= 3
  ) %>%
  pull(Genet)


within_pca <- lapply(
  genets_pca,
  run_within_genet_pca,
  geno = GT_num,
  genet_map = genet_map
)

names(within_pca) <- genets_pca


# Remove NULL results if any
within_pca <- within_pca[
  !vapply(
    within_pca,
    is.null,
    logical(1)
  )
]


# ================================================================
# 8. PCA SUMMARY
# ================================================================

within_pca_summary <- bind_rows(
  lapply(
    within_pca,
    function(x) {
      
      data.frame(
        Genet = x$Genet,
        N_samples = x$N_samples,
        N_variable_SNPs = x$N_variable_SNPs,
        PC1_percent = x$PC1_percent,
        PC2_percent = x$PC2_percent
      )
    }
  )
)

within_pca_summary <- within_pca_summary %>%
  arrange(
    as.numeric(
      sub("G", "", Genet)
    )
  )

print(
  within_pca_summary,
  n = Inf
)


# ================================================================
# 9. COMBINE PCA COORDINATES
# ================================================================

within_pca_scores <- bind_rows(
  lapply(
    within_pca,
    function(x) x$scores
  )
)


# ================================================================
# 10. OPTIONAL: ADD SAMPLE TYPE
#
# If genet_map already contains Sample_type, this will add it.
# Otherwise the PCA can be plotted without Sample_type.
# ================================================================

if ("Sample_type" %in% names(genet_map)) {
  
  within_pca_scores <- within_pca_scores %>%
    left_join(
      genet_map %>%
        select(
          Sample,
          Sample_type
        ) %>%
        distinct(),
      by = "Sample"
    )
}


# ================================================================
# 11. PLOT ALL WITHIN-GENET PCA RESULTS
# ================================================================

if ("Sample_type" %in% names(within_pca_scores)) {
  
  p_within_genet <- ggplot(
    within_pca_scores,
    aes(
      x = PC1,
      y = PC2,
      shape = Sample_type
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ Genet,
      scales = "free"
    ) +
    theme_classic(
      base_family = "serif"
    ) +
    labs(
      x = "PC1",
      y = "PC2",
      shape = "Sample type"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      )
    )
  
} else {
  
  p_within_genet <- ggplot(
    within_pca_scores,
    aes(
      x = PC1,
      y = PC2
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ Genet,
      scales = "free"
    ) +
    theme_classic(
      base_family = "serif"
    ) +
    labs(
      x = "PC1",
      y = "PC2"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      )
    )
}

p_within_genet


# ================================================================
# PART B
# EMPIRICAL ALLELE-BALANCE DISTRIBUTIONS
# ================================================================


# ================================================================
# 12. ALLELE BALANCE MATRIX
# ================================================================

Allele_balance <-
  vcf_fields$AD_ALT /
  (
    vcf_fields$AD_REF +
      vcf_fields$AD_ALT
  )


# ================================================================
# 13. EMPIRICAL AB DISTRIBUTION:
#     HETEROZYGOUS CALLS
# ================================================================

het_AB_all <- data.frame(
  GT = as.vector(vcf_fields$GT),
  DP = as.vector(vcf_fields$DP),
  GQ = as.vector(vcf_fields$GQ),
  AD_REF = as.vector(vcf_fields$AD_REF),
  AD_ALT = as.vector(vcf_fields$AD_ALT)
) %>%
  filter(
    GT %in% c("0/1", "1/0"),
    DP >= 10,
    GQ >= 30
  ) %>%
  mutate(
    AB =
      AD_ALT /
      (AD_REF + AD_ALT)
  )


summary(
  het_AB_all$AB
)

quantile(
  het_AB_all$AB,
  probs = c(
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99
  ),
  na.rm = TRUE
)


het_AB_summary <- het_AB_all %>%
  summarise(
    
    N = n(),
    
    Median_AB =
      median(
        AB,
        na.rm = TRUE
      ),
    
    Q25_AB =
      quantile(
        AB,
        0.25,
        na.rm = TRUE
      ),
    
    Q75_AB =
      quantile(
        AB,
        0.75,
        na.rm = TRUE
      ),
    
    Percent_020_080 =
      100 * mean(
        AB >= 0.20 &
          AB <= 0.80,
        na.rm = TRUE
      ),
    
    Percent_025_075 =
      100 * mean(
        AB >= 0.25 &
          AB <= 0.75,
        na.rm = TRUE
      ),
    
    Percent_both3 =
      100 * mean(
        AD_REF >= 3 &
          AD_ALT >= 3,
        na.rm = TRUE
      )
  )

het_AB_summary


# ================================================================
# 14. EMPIRICAL AB DISTRIBUTION:
#     HOMOZYGOUS CALLS
# ================================================================

hom_AB_all <- data.frame(
  GT = as.vector(vcf_fields$GT),
  DP = as.vector(vcf_fields$DP),
  GQ = as.vector(vcf_fields$GQ),
  AD_REF = as.vector(vcf_fields$AD_REF),
  AD_ALT = as.vector(vcf_fields$AD_ALT)
) %>%
  filter(
    GT %in% c("0/0", "1/1"),
    DP >= 10,
    GQ >= 30
  ) %>%
  mutate(
    AB =
      AD_ALT /
      (AD_REF + AD_ALT)
  )


hom_AB_summary <- hom_AB_all %>%
  group_by(GT) %>%
  summarise(
    
    N = n(),
    
    Median_AB =
      median(
        AB,
        na.rm = TRUE
      ),
    
    Q01 =
      quantile(
        AB,
        0.01,
        na.rm = TRUE
      ),
    
    Q05 =
      quantile(
        AB,
        0.05,
        na.rm = TRUE
      ),
    
    Q95 =
      quantile(
        AB,
        0.95,
        na.rm = TRUE
      ),
    
    Q99 =
      quantile(
        AB,
        0.99,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

hom_AB_summary


# ================================================================
# PART C
# WITHIN-GENET GENOTYPE DISCORDANCES
# ================================================================


# ================================================================
# 15. FUNCTION TO ANALYSE ONE GENET
#
# Final QC:
#
# DP >= 10
# GQ >= 30
#
# 0/0: AB <= 0.05
# 0/1: 0.25 <= AB <= 0.75
# 1/1: AB >= 0.95
# ================================================================

analyse_genet_global <- function(genet_name) {
  
  samples <- genet_map$Sample[
    genet_map$Genet == genet_name
  ]
  
  samples <- intersect(
    samples,
    colnames(vcf_fields$GT)
  )
  
  if (length(samples) < 3) {
    return(NULL)
  }
  
  gt <- vcf_fields$GT[
    ,
    samples,
    drop = FALSE
  ]
  
  gt_num <- matrix(
    NA_real_,
    nrow = nrow(gt),
    ncol = ncol(gt),
    dimnames = dimnames(gt)
  )
  
  gt_num[gt == "0/0"] <- 0
  
  gt_num[
    gt %in% c("0/1", "1/0")
  ] <- 1
  
  gt_num[gt == "1/1"] <- 2
  
  
  # Consensus genotype for each SNP
  consensus <- apply(
    gt_num,
    1,
    consensus_gt
  )
  
  
  consensus_text <- ifelse(
    consensus == 0,
    "0/0",
    ifelse(
      consensus == 1,
      "0/1",
      ifelse(
        consensus == 2,
        "1/1",
        NA_character_
      )
    )
  )
  
  
  res <- list()
  
  
  for (s in samples) {
    
    x <- gt_num[, s]
    
    idx <- which(
      !is.na(x) &
        !is.na(consensus) &
        x != consensus
    )
    
    if (length(idx) == 0) {
      next
    }
    
    
    gt_sample <-
      vcf_fields$GT[
        idx,
        s
      ]
    
    dp <-
      vcf_fields$DP[
        idx,
        s
      ]
    
    gq <-
      vcf_fields$GQ[
        idx,
        s
      ]
    
    ad_ref <-
      vcf_fields$AD_REF[
        idx,
        s
      ]
    
    ad_alt <-
      vcf_fields$AD_ALT[
        idx,
        s
      ]
    
    
    ab <-
      ad_alt /
      (ad_ref + ad_alt)
    
    
    DP_GQ_OK <-
      !is.na(dp) &
      !is.na(gq) &
      dp >= 10 &
      gq >= 30
    
    
    Allelic_support_OK <- case_when(
      
      gt_sample == "0/0" ~
        !is.na(ab) &
        ab <= 0.05,
      
      gt_sample %in% c(
        "0/1",
        "1/0"
      ) ~
        !is.na(ab) &
        ab >= 0.25 &
        ab <= 0.75,
      
      gt_sample == "1/1" ~
        !is.na(ab) &
        ab >= 0.95,
      
      TRUE ~ FALSE
    )
    
    
    res[[s]] <- data.frame(
      
      Genet = genet_name,
      
      N_samples_genet =
        length(samples),
      
      Sample = s,
      
      SNP_index = idx,
      
      SNP_ID =
        vcf$ID[idx],
      
      Chromosome =
        vcf$`#CHROM`[idx],
      
      Position =
        vcf$POS[idx],
      
      GT_consensus =
        consensus_text[idx],
      
      GT_sample =
        gt_sample,
      
      DP = dp,
      
      GQ = gq,
      
      AD_REF = ad_ref,
      
      AD_ALT = ad_alt,
      
      Allele_balance = ab,
      
      DP_GQ_OK =
        DP_GQ_OK,
      
      Allelic_support_OK =
        Allelic_support_OK,
      
      Final_QC =
        DP_GQ_OK &
        Allelic_support_OK,
      
      stringsAsFactors = FALSE
    )
  }
  
  
  bind_rows(res)
}


# ================================================================
# 16. GENETS INCLUDED IN CONSENSUS ANALYSIS
#
# Only genets with >= 3 samples.
#
# G8, G17, G23, G24 and G25 have only 2 samples
# and are therefore excluded.
# ================================================================

genets_to_analyse <- genet_map %>%
  count(
    Genet,
    name = "N_samples"
  ) %>%
  filter(
    N_samples >= 3
  ) %>%
  pull(Genet)

length(genets_to_analyse)

genets_to_analyse


# ================================================================
# 17. RUN ANALYSIS ACROSS ALL ELIGIBLE GENETS
# ================================================================

all_genet_list <- lapply(
  genets_to_analyse,
  analyse_genet_global
)

names(all_genet_list) <-
  genets_to_analyse


all_genet_discordances <- bind_rows(
  all_genet_list
)


dim(
  all_genet_discordances
)


# ================================================================
# 18. SUMMARY BY GENET
# ================================================================

final_genet_table <- all_genet_discordances %>%
  group_by(
    Genet,
    N_samples_genet
  ) %>%
  summarise(
    
    Raw_discordances =
      n(),
    
    Raw_loci =
      n_distinct(SNP_ID),
    
    After_DP_GQ =
      sum(
        DP_GQ_OK,
        na.rm = TRUE
      ),
    
    Loci_after_DP_GQ =
      n_distinct(
        SNP_ID[DP_GQ_OK]
      ),
    
    After_final_QC =
      sum(
        Final_QC,
        na.rm = TRUE
      ),
    
    Final_loci =
      n_distinct(
        SNP_ID[Final_QC]
      ),
    
    Percent_final =
      round(
        100 *
          After_final_QC /
          After_DP_GQ,
        1
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    as.numeric(
      sub("G", "", Genet)
    )
  )


print(
  final_genet_table,
  n = Inf
)


# ================================================================
# 19. GLOBAL QC SUMMARY
# ================================================================

global_QC_summary <- all_genet_discordances %>%
  summarise(
    
    N_genets =
      n_distinct(Genet),
    
    Raw_discordances =
      n(),
    
    Raw_loci =
      n_distinct(SNP_ID),
    
    After_DP_GQ =
      sum(
        DP_GQ_OK,
        na.rm = TRUE
      ),
    
    After_final_QC =
      sum(
        Final_QC,
        na.rm = TRUE
      ),
    
    Final_loci_global =
      n_distinct(
        SNP_ID[Final_QC]
      ),
    
    Percent_final =
      round(
        100 *
          After_final_QC /
          After_DP_GQ,
        1
      )
  )


global_QC_summary


# ================================================================
# PART D
# RECURRENCE OF HIGH-QUALITY DISCORDANCES
# ================================================================


# ================================================================
# 20. NUMBER OF SAMPLES SHARING EACH DISCORDANT GENOTYPE
# ================================================================

recurrent_final <- all_genet_discordances %>%
  filter(
    Final_QC
  ) %>%
  group_by(
    Genet,
    SNP_ID,
    Chromosome,
    Position,
    GT_consensus,
    GT_sample
  ) %>%
  summarise(
    
    N_samples_with_GT =
      n_distinct(Sample),
    
    Samples =
      paste(
        sort(
          unique(Sample)
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
  )


# ================================================================
# 21. RECURRENCE SUMMARY BY GENET
# ================================================================

recurrence_summary <- recurrent_final %>%
  group_by(Genet) %>%
  summarise(
    
    High_quality_loci =
      n_distinct(SNP_ID),
    
    Recurrent_2 =
      n_distinct(
        SNP_ID[
          N_samples_with_GT >= 2
        ]
      ),
    
    Recurrent_5 =
      n_distinct(
        SNP_ID[
          N_samples_with_GT >= 5
        ]
      ),
    
    Recurrent_10 =
      n_distinct(
        SNP_ID[
          N_samples_with_GT >= 10
        ]
      ),
    
    Max_samples =
      max(
        N_samples_with_GT,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    as.numeric(
      sub("G", "", Genet)
    )
  )


print(
  recurrence_summary,
  n = Inf
)


# ================================================================
# 22. GLOBAL RECURRENCE SUMMARY
#
# IMPORTANT:
# count unique GENET x SNP combinations
# ================================================================

global_recurrence <- recurrent_final %>%
  summarise(
    
    Total_high_quality_loci =
      n_distinct(
        paste(
          Genet,
          SNP_ID
        )
      ),
    
    Recurrent_2 =
      n_distinct(
        paste(
          Genet[
            N_samples_with_GT >= 2
          ],
          SNP_ID[
            N_samples_with_GT >= 2
          ]
        )
      ),
    
    Recurrent_5 =
      n_distinct(
        paste(
          Genet[
            N_samples_with_GT >= 5
          ],
          SNP_ID[
            N_samples_with_GT >= 5
          ]
        )
      ),
    
    Recurrent_10 =
      n_distinct(
        paste(
          Genet[
            N_samples_with_GT >= 10
          ],
          SNP_ID[
            N_samples_with_GT >= 10
          ]
        )
      )
  )


global_recurrence


# ================================================================
# PART E
# LOCAL RECURRENT GENOTYPIC SIGNATURES
# ================================================================


# ================================================================
# 23. LOCI RECURRENT IN >= 5 SAMPLES
# ================================================================

recurrent_5 <- recurrent_final %>%
  filter(
    N_samples_with_GT >= 5
  ) %>%
  arrange(
    Genet,
    Chromosome,
    Position
  )


nrow(recurrent_5)


# ================================================================
# 24. GROUP RECURRENT SNPs INTO LOCAL CLUSTERS
#
# Operational clustering threshold:
# consecutive SNPs <= 100 bp apart
# ================================================================

recurrent_5_clusters <- recurrent_5 %>%
  group_by(
    Genet,
    Chromosome
  ) %>%
  arrange(
    Position,
    .by_group = TRUE
  ) %>%
  mutate(
    
    Gap =
      Position -
      lag(Position),
    
    New_cluster =
      ifelse(
        is.na(Gap) |
          Gap > 100,
        1L,
        0L
      ),
    
    Cluster_number =
      cumsum(
        New_cluster
      )
  ) %>%
  ungroup()


# ================================================================
# 25. SUMMARIZE PHYSICAL CLUSTERS
# ================================================================

cluster_summary <- recurrent_5_clusters %>%
  group_by(
    Genet,
    Chromosome,
    Cluster_number
  ) %>%
  summarise(
    
    Start =
      min(Position),
    
    End =
      max(Position),
    
    Span_bp =
      End - Start,
    
    N_SNPs =
      n_distinct(SNP_ID),
    
    Max_samples =
      max(
        N_samples_with_GT
      ),
    
    SNPs =
      paste(
        SNP_ID,
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    Genet,
    Chromosome,
    Start
  )


# ================================================================
# 26. MULTI-SNP VS ISOLATED RECURRENT LOCI
# ================================================================

multilocus_signatures <- cluster_summary %>%
  filter(
    N_SNPs >= 2
  )


isolated_recurrent_loci <- cluster_summary %>%
  filter(
    N_SNPs == 1
  )


nrow(cluster_summary)

nrow(multilocus_signatures)

nrow(isolated_recurrent_loci)

sum(
  multilocus_signatures$N_SNPs
)


# ================================================================
# 27. VERIFY CO-OCCURRENCE IN THE SAME SAMPLES
#
# A physical cluster is considered a shared local genotypic
# signature only if the same samples carry the discordant genotype
# across all SNPs of the cluster.
# ================================================================

cluster_sample_overlap <- recurrent_5_clusters %>%
  group_by(
    Genet,
    Chromosome,
    Cluster_number
  ) %>%
  filter(
    n_distinct(SNP_ID) >= 2
  ) %>%
  summarise(
    
    N_SNPs =
      n_distinct(SNP_ID),
    
    Start =
      min(Position),
    
    End =
      max(Position),
    
    Span_bp =
      End - Start,
    
    N_shared_samples = {
      
      sample_lists <- strsplit(
        Samples,
        "; ",
        fixed = TRUE
      )
      
      length(
        Reduce(
          intersect,
          sample_lists
        )
      )
    },
    
    Shared_samples = {
      
      sample_lists <- strsplit(
        Samples,
        "; ",
        fixed = TRUE
      )
      
      x <- Reduce(
        intersect,
        sample_lists
      )
      
      paste(
        x,
        collapse = "; "
      )
    },
    
    .groups = "drop"
  )


# ================================================================
# 28. SUMMARY OF SHARED LOCAL SIGNATURES
# ================================================================

local_signature_summary <- cluster_sample_overlap %>%
  summarise(
    
    N_clusters =
      n(),
    
    N_with_shared_samples =
      sum(
        N_shared_samples >= 1
      ),
    
    N_shared_2 =
      sum(
        N_shared_samples >= 2
      ),
    
    N_shared_5 =
      sum(
        N_shared_samples >= 5
      ),
    
    Max_shared =
      max(
        N_shared_samples
      )
  )


local_signature_summary


table(
  cluster_sample_overlap$N_shared_samples
)


# ================================================================
# 29. EXAMPLE: G6 CHROMOSOME 15
# ================================================================

cluster_sample_overlap %>%
  filter(
    Genet == "G6",
    Chromosome == "Potre.1MX.Chr15"
  )


# ================================================================
# PART F
# FINAL SUPPLEMENTARY TABLE
# ================================================================


# ================================================================
# 30. COMBINE QC + RECURRENCE
# ================================================================

supplementary_within_genet <- final_genet_table %>%
  select(
    Genet,
    N_samples = N_samples_genet,
    Raw_discordances,
    After_DP_GQ,
    After_final_QC
  ) %>%
  left_join(
    recurrence_summary %>%
      select(
        Genet,
        High_quality_loci,
        Recurrent_2,
        Recurrent_5,
        Recurrent_10,
        Max_samples
      ),
    by = "Genet"
  ) %>%
  arrange(
    as.numeric(
      sub("G", "", Genet)
    )
  )


print(
  supplementary_within_genet,
  n = Inf
)


# ================================================================
# PART G
# EXPORT RESULTS
# ================================================================

write_csv(
  within_pca_summary,
  "within_genet_PCA_summary.csv"
)

write_csv(
  within_pca_scores,
  "within_genet_PCA_coordinates.csv"
)

write_csv(
  global_QC_summary,
  "within_genet_global_QC_summary.csv"
)

write_csv(
  supplementary_within_genet,
  "within_genet_summary_by_genet.csv"
)

write_csv(
  recurrent_final,
  "within_genet_recurrent_loci.csv"
)

write_csv(
  cluster_summary,
  "within_genet_local_clusters.csv"
)

write_csv(
  cluster_sample_overlap,
  "within_genet_multilocus_signatures.csv"
)

write_csv(
  local_signature_summary,
  "within_genet_local_signature_summary.csv"
)


# ================================================================
# 31. SAVE ALL IMPORTANT R OBJECTS
# ================================================================

save(
  within_pca,
  within_pca_summary,
  within_pca_scores,
  all_genet_discordances,
  final_genet_table,
  global_QC_summary,
  recurrent_final,
  recurrence_summary,
  global_recurrence,
  recurrent_5,
  cluster_summary,
  multilocus_signatures,
  cluster_sample_overlap,
  local_signature_summary,
  supplementary_within_genet,
  file = "within_genet_complete_analysis.RData"
)


# ================================================================
# END
# ================================================================