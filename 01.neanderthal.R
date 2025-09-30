.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))
library(data.table)

# df <- fread("/scratch/richards/chen-yang.su/data/30.sgdp_papuans/untar/CEU.chr12.ND_match")
# df1 <- fread("/scratch/richards/chen-yang.su/data/30.sgdp_papuans/untar/IBS.chr12.ND_match")

eur <- NULL
for (chr in 1:22) {
  d  <- fread(paste0("/scratch/richards/chen-yang.su/data/30.sgdp_papuans/untar/GBR.chr", chr, ".ND_match"))
  eur <- rbind(eur, d)
}

eas <- NULL
for (chr in 1:22) {
  d  <- fread(paste0("/scratch/richards/chen-yang.su/data/30.sgdp_papuans/untar/JPT.chr", chr, ".ND_match"))
  eas <- rbind(eas, d)
}

eas1 <- NULL
for (chr in 1:22) {
  d  <- fread(paste0("/scratch/richards/chen-yang.su/data/30.sgdp_papuans/untar/KHV.chr", chr, ".ND_match"))
  eas1 <- rbind(eas, d)
}


setwd('/scratch/richards/chen-yang.su/01.pQTL/')
ea <- fread("01.aric/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
ed <- fread("01.decode/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
ef <- fread("01.fenland/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
eu <- fread("01.ukbppp/03.cis_pQTL/all_strict_v2g_cis_res.tsv") %>% rename("uniprot" = "uniprotid")
aa <- fread("01.aric_afr_5e-8/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
au <- fread("01.ukbppp_afr_5e-8/03.cis_pQTL/all_strict_v2g_cis_res.tsv") %>% rename("uniprot" = "uniprotid")
ek <- fread("01.kyoto_5e-8/03.cis_pQTL/all_strict_v2g_cis_res.tsv") 
eukb <- fread("/scratch/richards/chen-yang.su/data/01.1.format_ukbppp_eas/01.ukbppp_eas_5e-8/03.cis_pQTL/all_strict_v2g_cis_res.tsv") 




# Check overlap for EUR and EAS

# 1. Number and percentage of pQTLs that overlap 
overlap1 <- eur[eur$ID %in% ea$rs_id, ]  # 20
overlap2 <- eur[eur$ID %in% ed$rs_id, ]  # 12
overlap3 <- eur[eur$ID %in% ef$rs_id, ]  # 21
overlap4 <- eur[eur$ID %in% eu$rs_id, ]  # 10

overlap5_JPT <- eas[eas$ID %in% ek$rs_id, ]  # Kyoto EAS has 8 hits with JPT file
overlap5_KHV <- eas1[eas1$ID %in% ek$rs_id, ]  # Kyoto EAS has 8 hits with KHV file

overlap6_JPT <- eas[eas$ID %in% eukb$rs_id, ]  # UKB EAS has none with JPT file
overlap6_KHV <- eas1[eas1$ID %in% eukb$rs_id, ]  # UKB EAS has none with KHV file


# Summary of checking on alleles across 4 EUR cohorts: only 2 rsid removed from decode for genes SDSL and BCAT2 which are not relevant
###################################
# Updated overlap but match on rsid + alleles also (ensure REF/ALT matches pQTL effect/other or vice versa

# Ensure consistent column names for joining (optional but can make it clearer)
# For example, if eur$ID and ea$rs_id are the common identifiers
eur_clean <- eur %>% rename(rsid = ID)


ea_clean <- ea %>% rename(rsid = rs_id)
# Perform the join, then filter for allele matches
final_overlap1 <- inner_join(eur_clean, ea_clean, by = "rsid") %>%
  filter(
    (REF == effect_allele & ALT == other_allele) |
      (REF == other_allele & ALT == effect_allele)
  )

ed_clean <- ed %>% rename(rsid = rs_id)
# Perform the join, then filter for allele matches
final_overlap2 <- inner_join(eur_clean, ed_clean, by = "rsid") %>%
  filter(
    (REF == effect_allele & ALT == other_allele) |
      (REF == other_allele & ALT == effect_allele)
  )

ef_clean <- ef %>% rename(rsid = rs_id)
# Perform the join, then filter for allele matches
final_overlap3 <- inner_join(eur_clean, ef_clean, by = "rsid") %>%
  filter(
    (REF == effect_allele & ALT == other_allele) |
      (REF == other_allele & ALT == effect_allele)
  )

eu_clean <- eu %>% rename(rsid = rs_id)
# Perform the join, then filter for allele matches
final_overlap4 <- inner_join(eur_clean, eu_clean, by = "rsid") %>%
  filter(
    (REF == effect_allele & ALT == other_allele) |
      (REF == other_allele & ALT == effect_allele)
  )


# Any rsid got removed?
overlap1[!overlap1$ID %in% final_overlap1$rsid,]
overlap2[!overlap2$ID %in% final_overlap2$rsid,]
# > overlap2[!overlap2$ID %in% final_overlap2$rsid,]
# CHROM       POS         ID    REF    ALT SEGMENT ALLELE  SCORE   NMATCH   DMATCH
# <int>     <int>     <char> <char> <char>   <int>  <int>  <int>   <char>   <char>
#   1:    12 113835773 rs71467906      G      C       3      1 684305 mismatch mismatch
# 2:    19  49304059  rs4801775      T      C      14      1 237577  notcomp  notcomp
ed[ed$rs_id == "rs71467906",]
ed[ed$rs_id == "rs4801775",]


overlap3[!overlap3$ID %in% final_overlap3$rsid,]
overlap4[!overlap4$ID %in% final_overlap4$rsid,]

###################################


# left merge strict v2g gene
overlap1 <- merge(overlap1, ea[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")
overlap2 <- merge(overlap2, ed[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")
overlap3 <- merge(overlap3, ef[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")
overlap4 <- merge(overlap4, eu[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")

overlap5_JPT <- merge(overlap5_JPT, ek[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")
overlap5_KHV <- merge(overlap5_KHV, ek[,c("rs_id", "gene")], by.x="ID", by.y="rs_id", all.x=T) %>% rename("instrumented_gene" = "gene")


# summary ----

# Generate cross-tabulation

get_summary <-function(overlap) {
  print(nrow(overlap))  # number of strict v2g cis-pQTLs that have an rsid matching in the Neanderthal file
  print(overlap)  # print these overlapping variants
  print(overlap[overlap$NMATCH=="match"])  # print only the ones with NMATCH=match
  combo_table <- table(overlap$NMATCH, overlap$DMATCH)
  combo_df <- as.data.frame(combo_table) 
  colnames(combo_df) <- c("Neanderthal", "Denisovan", "Count")
  print(combo_df)
  print("*********************************")
}

get_summary(overlap1)
get_summary(overlap2)
get_summary(overlap3)
get_summary(overlap4)


table(c(overlap1[overlap1$NMATCH=="match"]$instrumented_gene, 
      overlap2[overlap2$NMATCH=="match"]$instrumented_gene,
      overlap3[overlap3$NMATCH=="match"]$instrumented_gene,
      overlap4[overlap4$NMATCH=="match"]$instrumented_gene))


# Overlap Fenland and UKB-PPP
table(c(overlap3[overlap3$NMATCH=="match"]$instrumented_gene,
        overlap4[overlap4$NMATCH=="match"]$instrumented_gene))
# ADAMTS13   ADIPOQ     BCL2     BMP1    CCL25   COL9A1     FRZB  HTATIP2   IL11RA   LGALS3     LHPP  LRRFIP1    MMP12    NELL1     NTN4     OAS1 
# 1        1        1        1        1        1        1        1        1        1        1        1        1        2        1        1 
# PLXDC2 TNFRSF1B 
# 1        1

# NELL1 appears twice
# Mention OAS1 positive control




overlap3[overlap3$NMATCH=="match"]
overlap4[overlap4$NMATCH=="match"]

stop()
###########





print(length(ea$rs_id))


ea[ea$rs_id %in% eur$ID,]
ed[ed$rs_id %in% eur$ID,]
ef[ef$rs_id %in% eur$ID,]
eu[eu$rs_id %in% eur$ID,]

eukb[eukb$rs_id %in% eas$ID,]
eukb[eukb$rs_id %in% eas1$ID,]










length(unique(df$SEGMENT))

# This SNP (rs4767027) at position 113,359,157 on chromosome 12:
#   
# Has a reference allele T and alternate allele C.
# 
# Belongs to inferred segment 2, with the introgressed allele predicted to be the reference (T).
# 
# The segment has a high score of 696811, indicating strong evidence of archaic introgression.
# 
# This SNP matches Neanderthal but does not match Denisovan.
# 
# This implies that this particular allele likely originated from a Neanderthal introgression event in the CEU (European) population.





# Proteomics overlap ----

# UKB EUR proteomics pQTL


