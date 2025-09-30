library(data.table)
mvpeur <- fread("/scratch/richards/chen-yang.su/data/14.MVP/gwas_paths/EUR_MVP_table.txt")

mvpafr <- fread("/scratch/richards/chen-yang.su/data/14.MVP/gwas_paths/AFR_MVP_table.txt")
mvpafr[mvpafr$Trait == "CircCAD",]$Description
# [1] "Coronary Artery / Coronary Heart Disease (including heart attack, angina)"



# BMP1 Fenland pQTL rs149714437
# 8:22198746 (GRCh38)
# 8:22056259 (GRCh37)
mvpafrgwas <- fread("/scratch/richards/public/MVP/GIA/MVP_R4.1000G_AGR.GIA.Surveys_Binary_batch1/MVP_R4.1000G_AGR.CircCAD.AFR.GIA.dbGaP.txt.gz")

mvpafrgwas_subset <- mvpafrgwas[mvpafrgwas$chrom == 8 & 
                           mvpafrgwas$pos >= 22198746 - 10000000 &  # pad 10Mb upstream and downstream of pQTL
                           mvpafrgwas$pos <= 22198746 + 10000000 
]

mvpafrgwas_subset <- mvpafrgwas_subset[order(mvpafrgwas_subset$pos),]             
write_tsv(mvpafrgwas_subset, "/scratch/richards/chen-yang.su/24.atlas/03.haplotype/CircCAD_regionGRCh38_MVPAFR_locuszoom.tsv")



mvp <- fread("/scratch/richards/chen-yang.su/13.pqtl_mvp/02.summarize_all/ST_mr_coloc_res/EUR_mr_coloc.tsv")

bmp1 <- mvp[mvp$protein == "BMP1",]
bmp1

bmp1$Description
# [1] "Coronary Artery / Coronary Heart Disease (including heart attack, angina)"
# [2] "Ischemic Heart Disease"                                                   
# [3] "Coronary atherosclerosis"                                                 
# [4] "Other chronic ischemic heart disease, unspecified"