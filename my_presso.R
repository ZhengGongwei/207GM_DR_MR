my_presso <- function(RISKlist,DRlist) {
    for(i in RISKlist){
    file1 <- paste0(i,"_1kr0.1.txt")
    exposure <- read.table(file1,sep="\t",header=T)
    exposure$id.exposure <- i
    for(j in DRlist){
        outcome <- extract_outcome_data(snps = exposure$SNP,outcomes = j)
        dat <- harmonise_data(exposure_dat = exposure,outcome_dat = outcome)
        dat2 <- dat[dat$mr_keep,]
        if(length(dat2$SNP) == 1){
            method_list1 <- c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
            res <- mr(dat, method_list = method_list1)
            het <- mr_heterogeneity(dat,method_list = method_list1[1])
            plt <- mr_pleiotropy_test(dat)
        } else if (length(dat2$SNP) > 1 & length(dat2$SNP) <= 3){
            method_list2 <- c("mr_ivw_fe","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
            res <- mr(dat, method_list = method_list2)
            het <- mr_heterogeneity(dat,method_list = method_list2[1])
            plt <- mr_pleiotropy_test(dat)
        } else if (length(dat2$SNP) > 3){
            method_list3 <- c("mr_ivw_mre","mr_ivw","mr_egger_regression","mr_two_sample_ml","mr_weighted_median","mr_weighted_mode")
            res <- mr(dat, method_list = method_list3)
            het <- mr_heterogeneity(dat,method_list = method_list3[1])
            plt <- mr_pleiotropy_test(dat)
        }

        if(dim(het)[2] == 8){
            het <- het[,5:8]
            # colnames(het) <- lapply(colnames(het), FUN=function(x) paste0("het_",x))
            aa <- merge(res[1,],het,by.x="method",by.y="method",all.x=TRUE)
        }else{
            het <- data.frame(res$method[1],NA,NA,NA)
            colnames(het) <- c("method","Q","Q_df","Q_pval")
            aa <- merge(res[1,],het,by.x="method",by.y="method",all.x=TRUE)
        }


        if(dim(plt)[2] == 7){
            plt <- plt[,c(5,7)]
            plt$method <- het$method
            colnames(plt) <- lapply(colnames(plt), FUN=function(x) paste0("plt_",x))
            aa <- merge(aa,plt,by.x="method",by.y="plt_method",all.x=TRUE)
        }else{
            plt <- data.frame(res$method[1],NA,NA)
            colnames(plt) <- c("het_method","plt_egger_intercept","plt_pval")
            aa <- merge(aa,plt,by.x="method",by.y="het_method",all.x=TRUE)
        }
            aa <- aa[,-c(4,5)]
    
            res_pre = tryCatch({res_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdOutcome = "se.outcome", SdExposure = "se.exposure",OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat2,NbDistribution = 1000,  SignifThreshold = 0.05)
                             data.frame("mr_presso_beta" = res_presso$`Main MR results`$`Causal Estimate`[1],
                                        "mr_presso_sd" = res_presso$`Main MR results`$`Sd`[1],
                                        "mr_presso_P-value" = res_presso$`Main MR results`$`P-value`[1],
                                        "GlobalTest_Pval" = res_presso$`MR-PRESSO results`$`Global Test`$Pvalue,
                                        "Outlier_nsnps" = length(res_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                                        "Outlier_corrected_beta" = res_presso$`Main MR results`$`Causal Estimate`[2],
                                        "Outlier_corrected_sd" = res_presso$`Main MR results`$`Sd`[2],
                                        "Outlier_corrected_P-value" = res_presso$`Main MR results`$`P-value`[2],
                                        "Distortion_P_value" = ifelse(is.null(res_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue), NA, res_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue))
                            }, warning = function(w) {
                                data.frame("mr_presso_beta" = res_presso$`Main MR results`$`Causal Estimate`[1],
                                        "mr_presso_sd" = res_presso$`Main MR results`$`Sd`[1],
                                        "mr_presso_P-value" = res_presso$`Main MR results`$`P-value`[1],
                                        "GlobalTest_Pval" = res_presso$`MR-PRESSO results`$`Global Test`$Pvalue,
                                        "Outlier_nsnps" = length(res_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                                        "Outlier_corrected_beta" = res_presso$`Main MR results`$`Causal Estimate`[2],
                                        "Outlier_corrected_sd" = res_presso$`Main MR results`$`Sd`[2],
                                        "Outlier_corrected_P-value" = res_presso$`Main MR results`$`P-value`[2],
                                        "Distortion_P_value" = ifelse(is.null(res_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue), NA, res_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue))
                            }, error = function(e) {
                                data.frame("mr_presso_beta" = NA,"mr_presso_sd" = NA,"mr_presso_P-value" = NA,"GlobalTest_Pval" = NA,"Outlier_nsnps" = NA,"Outlier_corrected_beta" = NA,
                                "Outlier_corrected_sd" = NA,"Outlier_corrected_P-value" = NA,"Distortion_P_value" = NA
                                )
                            }, finally = {
                            })
        aa <- cbind(aa,res_pre)
        #data2 <- rbind(data2,aa)
        fn <- paste0("C:/Users/Lenovo/OneDrive/MR/tmp/",i,"_",j,".txt")
        write.table(aa, file=fn, row=FALSE, col=TRUE, qu=FALSE, sep="\t")
        }
    }
}