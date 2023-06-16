#  pos /share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps/
#  所有taxa分类的菌

# env MR
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
setwd("/share2/pub/GM_GWAS_LD_clumped_snps/")

for(dir in c("/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_genus_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_family_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_order_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_phyllum_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_cluster_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_species_01022021_v1/"))
{

  expo_gwas_list <- list.files(path=dir,pattern = "*.tsv$")
  for(i in expo_gwas_list){
    indexs <- strsplit(i, ".", fixed= T)
    expo_name = indexs[[1]][length(indexs[[1]])-1]
    expogwas <- paste(dir,i,sep='')
    expo_dat <- read_exposure_data(
      filename = expogwas,
      sep = "\t",
      snp_col = "id",
      beta_col = "beta",
      se_col = "SE",
      effect_allele_col = "alt",#alt
      other_allele_col = "ref",#ref
      eaf_col = "AF_Allele2",
      pval_col = "pval"
    )
    expo_dat <- expo_dat[which(expo_dat$pval.exposure<1e-5),]
    # remove genus's LD
    b <- ld_clump(
    dplyr::tibble(rsid=expo_dat$SNP, pval=expo_dat$pval.exposure, id=expo_dat$id.exposure),
    plink_bin = "/share2/pub/GM_GWAS_LD_clumped_snps/plink",
    bfile = "/share2/pub/GM_GWAS_LD_clumped_snps/EUR_ref/EUR",
    clump_kb = 1000,clump_r2 = 0.1
    )
    expo_dat <- expo_dat[which(expo_dat$SNP %in% b$rsid),]
    LD_clumped <- paste0("/share2/pub/GM_GWAS_LD_clumped_snps/",expo_name,"_1kr0.1.txt")
    write.table(expo_dat,file=LD_clumped)
    }
}
