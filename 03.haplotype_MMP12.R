.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))
library(data.table)
library(tidyverse)
library(dplyr)

fenland <- fread("/scratch/richards/chen-yang.su/01.pQTL/01.ukbppp/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
fenland[fenland$gene == "OAS1",]$protein  # OAS1.10361_25 
fenland[fenland$gene == "BMP1",]$protein  # BMP1.3348_49 
fenland[fenland$gene == "MMP12",]  
pqtl_chr <- fenland[fenland$gene == "MMP12",]$chromosome
pqtl_pos <- fenland[fenland$gene == "MMP12",]$base_pair_location
pqtl_rsid <- fenland[fenland$gene == "MMP12",]$rs_id

# UKBPPP Protein GWAS (hg38)
bmp1 <- fread("/scratch/richards/chen-yang.su/data/01.1.format_ukbppp/untar_dir/MMP12_P39900_OID21439_v1_Oncology.gz")
bmp1$Allele1 <- toupper(bmp1$Allele1)  # ALT
bmp1$Allele2 <- toupper(bmp1$Allele2)  # REF
bmp1_region <- bmp1[bmp1$chr == pqtl_chr & 
                      bmp1$pos >= pqtl_pos - 10000000 &  # pad 10Mb upstream and downstream of pQTL
                      bmp1$pos <= pqtl_pos + 10000000 
]
bmp1_region <- bmp1_region[order(bmp1_region$pos),]             
write_tsv(bmp1_region, "/scratch/richards/chen-yang.su/24.atlas/03.haplotype/mmp12_regionGRCh38_ukbppp.tsv")


# Haplotype Length? ----

###########
# Step 1: Create a region file
region <- data.frame(
  chr = paste0(pqtl_chr),
  start = pqtl_pos - 10000000, 
  end = pqtl_pos + 10000000    
)
write.table(region, file = "/scratch/richards/chen-yang.su/24.atlas/03.haplotype/MMP12_region_hg38_PLINK.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 2: Calculate LD with the target SNP
cmd <- paste0(
  "/scratch/richards/chen-yang.su/software/plink ",
  "--bfile /scratch/richards/satoshi.yoshiji/database/1KG/EUR/EUR ",
  "--chr ", pqtl_chr, " ",
  "--from-bp ", region$start, " ",
  "--to-bp ", region$end, " ",
  "--ld-snp ", pqtl_rsid, " ",
  "--r2 ",
  "--ld-window 99999 ",
  "--ld-window-kb 200 ",
  "--ld-window-r2 0 ",
  "--out /scratch/richards/chen-yang.su/24.atlas/03.haplotype/ld_output_mmp12"
)
system(cmd)   # runs it

# 🔹 Step 3: Filter SNPs with r² > 0.8
ld <- read.table("/scratch/richards/chen-yang.su/24.atlas/03.haplotype/ld_output_mmp12.ld", header=TRUE)
high_ld <- subset(ld, R2 >= 0.7)

# left_37 <- which(min(high_ld$BP_B))
# right_37 <- max(high_ld$BP_B)
left_rsid <- high_ld[which.min(high_ld$BP_B), "SNP_B"]
right_rsid <- high_ld[which.max(high_ld$BP_B), "SNP_B"]


# # Get GRCh38 positions
# library(httr)
# library(jsonlite)
# 
# # Function to get SNP position from Ensembl
# get_snp_position <- function(snp) {
#   url <- paste0("https://rest.ensembl.org/variation/human/", snp, "?content-type=application/json")
#   response <- fromJSON(content(GET(url), as = "text"))
#   
#   if (!is.null(response$mappings)) {
#     grch38_pos <- response$mappings[response$mappings$assembly_name == "GRCh38", ]
#     return(data.frame(SNP = snp, Chr = grch38_pos$seq_region_name, Position = grch38_pos$start))
#   } else {
#     return(data.frame(SNP = snp, Chr = NA, Position = NA))
#   }
# }
# 
# # List of SNPs
# snps <- c(left_rsid, right_rsid)
# 
# # Get positions for all SNPs
# snp_positions <- do.call(rbind, lapply(snps, get_snp_position))
# 
# # Display results
# print(snp_positions)  # GRCh38 positions


###########

# On locuszoom (GRCh37): chr8 22056435 to 22017699
# haplo_df <- data.frame(
#   chr = paste0("chr", pqtl_chr),
#   start = 22160186, # GRCh37: 22017699
#   end = 22198922    # GRCh37: 22056435
# )
haplo_df <- data.frame(
  chr = paste0("chr", pqtl_chr),
  start = min(snp_positions$Position), 
  end = max(snp_positions$Position)   
)
haplo_length <- haplo_df$end - haplo_df$start
print(haplo_length)
write.table(haplo_df, file = "/scratch/richards/chen-yang.su/24.atlas/03.haplotype/MMP12_region_hg38.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Recombination rate? ----
# Example:
system(paste0("~/projects/richards/share/bin/bedtools map -a /scratch/richards/chen-yang.su/data/28.recombRate/Hugo_region_hg38.bed -b /scratch/richards/chen-yang.su/data/28.recombRate/bedGraph/recomb1000GAvg.bedGraph -c 4 -o mean"))
system(paste0("~/projects/richards/share/bin/bedtools map -a /scratch/richards/chen-yang.su/24.atlas/03.haplotype/MMP12_region_hg38.bed -b /scratch/richards/chen-yang.su/data/28.recombRate/bedGraph/recomb1000GAvg.bedGraph -c 4 -o mean"))

# 0.4208720269



# P common ancestor ----

# Define the function
p_common_ancestor <- function(length_bp,
                              recomb_cM_Mb,
                              t_split_yrs = 550000,
                              archaic_sample_age_yrs = 50000,
                              gen_time_yrs = 29) {
  # Recombination rate: Morgans per bp per generation
  r <- (recomb_cM_Mb / 100) / 1e6
  
  # Total branch length in generations
  t_tot_gen <- (t_split_yrs + (t_split_yrs - archaic_sample_age_yrs)) / gen_time_yrs
  
  # Expected tract length
  L_exp <- 1 / (r * t_tot_gen)
  
  # Scaled length
  x <- length_bp / L_exp
  
  # Gamma survival function (k = 2)
  return(exp(-x) * (1 + x))
}


# # ---------- toy replication of Zeberg & Pääbo 2020 ----------
# toy <- data.frame(
#   region_id = c("core 49.4 kb", "full 333.8 kb"),
#   length_bp = c(49400, 333800),
#   recomb_cM_Mb = c(0.53, 0.53)
# )
# 
# toy <- toy %>%
#   rowwise() %>%
#   mutate(P_common_ancestor = p_common_ancestor(length_bp, recomb_cM_Mb)) %>%
#   ungroup()
# 
# print(toy)


# Test yourself ----
toy2 <- data.frame(
  region_id = c("MMP12 UKBPPP haplotype"),
  length_bp = c(haplo_length),
  recomb_cM_Mb = c(0.4208720269)
)

toy2 <- toy2 %>%
  rowwise() %>%
  mutate(P_common_ancestor = p_common_ancestor(length_bp, recomb_cM_Mb)) %>%
  ungroup()

print(toy2)


# Use leftmost and rightmost variants ====
# > print(toy2)
# # A tibble: 1 × 4
# region_id              length_bp recomb_cM_Mb P_common_ancestor
# <chr>                      <dbl>        <dbl>             <dbl>
#   1 BMP1 Fenland haplotype     38736         3.80          3.62e-22

# Use R2 >= 0.8 ====
# # A tibble: 1 × 4
# region_id              length_bp recomb_cM_Mb P_common_ancestor
# <chr>                      <int>        <dbl>             <dbl>
#   1 BMP1 Fenland haplotype      5154         4.30           0.00296


# Use R2 >= 0.7 ====
# > print(toy2)
# # A tibble: 1 × 4
# region_id              length_bp recomb_cM_Mb P_common_ancestor
# <chr>                      <int>        <dbl>             <dbl>
#   1 BMP1 Fenland haplotype       176         6.15             0.941





# 3. Extract MR results ----

# https://broad.io/protein_mr_atlas
# BMP1 has hits for CAD (Aragam 2022) only


# 4.	Fine-mapping; Disease: CAD – check EUR and AFR (using locuszoom for fine-mapping, single-ancestry [Hugo’s Nat Gen OAS1 paper])
# Prep for locuszoom

# BMP1 Fenland pQTL rs149714437
# 8:22198746 (GRCh38)
# 8:22056259 (GRCh37)
cad <- fread("/scratch/richards/chen-yang.su/data/04.prep_MR_outcome/gwas_sum_stats/3_CAD_aragam_NatGen/gwas_catalog/CAD.GCST90132314.catalog.tsv.gz")
cad_subset <- cad[cad$chromosome == pqtl_chr & 
                    cad$base_pair_location >= 22198746 - 10000000 &  # pad 10Mb upstream and downstream of pQTL
                    cad$base_pair_location <= 22198746 + 10000000 
                  ]
cad_subset <- cad_subset[order(cad_subset$base_pair_location),]             
write_tsv(cad_subset, "/scratch/richards/chen-yang.su/24.atlas/03.haplotype/cad_regionGRCh38_EUR_locuszoom.tsv")







