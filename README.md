# neanderthal-introgression

# 01.neanderthal.R
Neanderthal/Denisovan overlap with cis-pQTLs (R script)

This R script aggregates archaic-introgression SNP match files (per chromosome) for several populations, then checks how many cis-pQTL lead variants (from multiple proteomics cohorts) overlap those archaic segments.
It reports:

Counts/lists of overlapping rsIDs for EUR (GBR reference) and EAS (JPT & KHV references).

An allele-aware overlap for EUR (ensuring REF/ALT match the pQTL effect_allele/other_allele).

Cross-tabs of whether each overlapping variant matches Neanderthal and/or Denisovan calls.

Quick gene-level tallies of instrumented genes among the overlaps (e.g., flags like OAS1).

No files are written—the results print to the console as data tables and counts.

# 02.p_common_ancestor_example.R
Probability of Common Ancestor from Tract Lengths (R)

## What this script does

Implements a simple probability model to estimate whether an observed genomic tract shares a common ancestor predating the human–Neanderthal/Denisovan split (i.e., not introgressed) given its length and local recombination rate. It reproduces toy values from Zeberg & Pääbo (2020) and lets you test alternative parameters.

## Method (brief)

For a tract of length L (bp) and recombination r (Morgans/bp/generation), with total branch length t_tot (generations),

Expected tract length: L_exp = 1 / (r * t_tot)

Scaled length: x = L / L_exp

Survival prob. (Gamma, k=2):
P_common_ancestor = exp(-x) * (1 + x)

Defaults:

Split time t_split_yrs = 550,000

Archaic sample age archaic_sample_age_yrs = 50,000

Generation time gen_time_yrs = 29
t_tot_gen = (t_split + (t_split − archaic_age)) / gen_time


# 03.haplotype_MMP12.R
MMP12 locus: region extraction, LD span, recomb rate & P(common ancestor)

## What this script does (short)

Loads UKB-PPP cis-pQTLs, grabs the MMP12 sentinel (chromosome, base_pair_location, rs_id).

Extracts a ±10 Mb region from UKB-PPP protein GWAS around the MMP12 pQTL and writes it to TSV.

Computes LD (EUR 1KG reference) around the pQTL via PLINK, keeps SNPs with R² ≥ 0.7, and derives the haplotype span (left/right rsIDs → GRCh38 positions; Ensembl REST call is provided but commented).

Writes the span to a BED file and queries a recombination-rate bedGraph with bedtools map to get the mean cM/Mb across the span.

Uses a small function p_common_ancestor() (Gamma survival, k=2) to estimate the probability the span predates human–Neanderthal split (i.e., not introgressed) given the span length and recombination rate.

Optionally subsets a CAD GWAS (Aragam 2022) ±10 Mb around the pQTL for LocusZoom.

# 03.haplotype.R
BMP1 (Fenland): region slice, LD span, recomb rate & P(common ancestor)

## What this script does (brief)

Loads Fenland cis-pQTLs, extracts BMP1 sentinel (chromosome, base_pair_location, rs_id).

Slices a ±10 Mb window from the Fenland protein GWAS (GRCh37/hg19) and writes TSV.

Computes LD around the pQTL using PLINK (1KG EUR), keeps SNPs with R² ≥ 0.7, and defines the haplotype span by the left/right edge SNPs.

Converts span edges to GRCh38 via Ensembl REST, writes a BED file.

Uses bedtools map with a recombination bedGraph to get mean cM/Mb over the span.

Evaluates p_common_ancestor() (Gamma survival, k=2) to estimate the probability the span predates the human–Neanderthal split (i.e., not introgressed).

Optionally extracts CAD GWAS (Aragam 2022) ±10 Mb for LocusZoom.

# 04.mvp_bmp1.R
MVP AFR slice around BMP1 pQTL & BMP1 MR/coloc lookup

## What this script does (brief)

Loads MVP EUR/AFR GWAS path tables.

Confirms the AFR phenotype description for CircCAD (coronary artery/heart disease umbrella).

Extracts an AFR MVP GWAS window (±10 Mb) around the BMP1 Fenland pQTL at 8:22198746 (GRCh38) and writes a LocusZoom-ready TSV.

Looks up BMP1 rows in a precomputed MR/coloc summary to list associated cardiovascular outcomes.


# 28.recombRate/01.convert.sh
Compute mean recombination rate for a region (hg38)

## What this does

Converts UCSC recombination bigWig tracks to bedGraph, then uses bedtools map to compute the mean local recombination rate (cM/Mb) for a given BED region. Includes a sanity-check region from Hugo’s paper (lifted to hg38).
