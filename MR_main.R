#pos /share2/pub/zhenggw/zhenggw/TwoSampleMR/eyedisease/NC/207GM-DM_RETINOPATHY/MRresult/
# env MR
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
setwd("/share2/pub/zhenggw/zhenggw/TwoSampleMR/eyedisease/NC/207GM-DM_RETINOPATHY/MRresult/")

for(dir in c("/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_genus_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_family_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_order_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_phyllum_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_cluster_01022021_v1/",
              "/share2/pub/lijj/lijj/TheDutchMicrobiomeProject/gwas_taxa_species_01022021_v1/"))
{

  expo_gwas_list <- list.files(path=dir,pattern = "*.tsv$")
  out_gwas_list <- c("finn-b-DM_RETINOPATHY","finn-b-H7_RETINOPATHYDIAB","finn-b-DM_RETINOPATHY_EXMORE")

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
    expo_dat <- clump_data(expo_dat,clump_kb = 1000,clump_r2 = 0.1,pop = "EUR")

    for(j in out_gwas_list){
      out_name <- gsub("_","-",j)

    res = tryCatch({
      outcome_dat<- extract_outcome_data(snps = expo_dat$SNP,outcomes = j)
      dat <- harmonise_data(exposure_dat = expo_dat,outcome_dat = outcome_dat)
      dat2 <- dat[dat$mr_keep,]
        if(length(dat2$SNP) == 1){
          res <- mr(dat, method_list = c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        } else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
          res <- mr(dat, method_list = c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        } else if (length(dat2$SNP) > 3){
          res <- mr(dat, method_list = c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        }
          }, warning = function(w){
      outcome_dat<- extract_outcome_data(snps = expo_dat$SNP,outcomes = j)
      dat <- harmonise_data(exposure_dat = expo_dat,outcome_dat = outcome_dat)
      dat2 <- dat[dat$mr_keep,]
        if(length(dat2$SNP) == 1){
          res <- mr(dat, method_list = c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        } else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
          res <- mr(dat, method_list = c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        } else if (length(dat2$SNP) > 3){
          res <- mr(dat, method_list = c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
        }   
          }, error = function(e){print("error exist!")},finally = {print("completed!")}
          )

      filename <- paste(expo_name,"_",out_name,"_res.txt",sep='')
      write.table(res, file = filename, sep = "\t",col.names = TRUE,row.names = FALSE)
      print(filename)
      }
    }
}
