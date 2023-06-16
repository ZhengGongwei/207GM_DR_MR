# prepare full gwas data and format (7*7)
# /share2/pub/zhenggw/zhenggw/MVMR
```shell
for file in *.tsv; do 
    awk '$2 !~ /^X/{ 
        FS="\t"; 
        OFS="\t"; 
        if(NR==1) { 
            for(i=1;i<=NF;i++) { 
                if($i=="id") {id=i} 
                if($i=="beta") {beta=i} 
                if($i=="SE") {SE=i} 
                if($i=="AF_Allele2") {AF_Allele2=i} 
                if($i=="alt") {alt=i} 
                if($i=="ref") {ref=i} 
                if($i=="pval") {pval=i} 
                if($i=="num") {num=i} 
            } 
            print "SNP", "beta", "se", "eaf", "effect_allele", 
                "other_allele", "pval", "samplesize" 
        } 
        else { 
            if($id=="" || $beta=="" || $SE=="" || $AF_Allele2=="" || $alt=="" || $ref=="" || $pval=="" || $num=="") { 
                print "Error: Missing column(s) in line " NR > "/dev/stderr"; 
                exit 1 
            } 
            else { 
                print $id, $beta, $SE, $AF_Allele2, $alt, $ref, $pval, $num 
            } 
        } 
    }' "$file" > "${file%.tsv}.txt"
done

for file in *.txt
do
    sed -i '1s/$/\tPhenotype/' $file    # 在第一行末尾添加列名
    sed -i "2,\$s/$/\t${file%.*}/" $file    # 在第二行到最后一行末尾添加对应文件的名称
done

cd /share2/pub/zhenggw/zhenggw/MVMR/Risk_data
for file in ieu-b-38.vcf.gz ieu-b-39.vcf.gz ieu-b-110.vcf.gz; do
    output_file="${file%%.*}.txt"
    zless "$file" | grep -v "^##" | grep -v "^X" | cut -f 1-5,10 | sed 's/:/\t/g' | awk '{print $3,$6,$7,$9,$5,$4,$8,$10}' | awk 'BEGIN {OFS="\t"} NR==1 {print "SNP\tbeta\tse\teaf\teffect_allele\tother_allele\tpval\tsamplesize"} NR!=1 {$7=10^(-$7); print}' > "$output_file"
done

awk 'NR==1 {print $0} NR>1{$8="440546"; print}' ieu-b-110.txt > output.txt

# IL16_v3.txt
less IL16_v3.txt | awk -v OFS='\t' '{print $12,$6,$7,$4,$2,$3,$8,$10}' | sed '1i SNP\tbeta\tse\teaf\teffect_allele\tother_allele\tpval\tsamplesize' | sed '2d'  > IL16.txt

# MAGIC1000G_2hGlu_EUR.tsv.gz MAGIC1000G_FG_EUR.tsv.gz MAGIC1000G_HbA1c_EUR.tsv.gz
for file in MAGIC1000G_2hGlu_EUR.tsv.gz MAGIC1000G_FG_EUR.tsv.gz MAGIC1000G_HbA1c_EUR.tsv.gz;do
	output_file1=`echo $file | awk -F "_" '{print $2}'`
	output_file="$output_file1.txt"
	zless "$file" | awk '$2 !~ /^X/{print}' | awk -v OFS="\t" '{print $1,$7,$8,$6,$4,$5,$9,$10}' | sed '1i SNP\tbeta\tse\teaf\teffect_allele\tother_allele\tpval\tsamplesize' | sed '2d' > "$output_file"
done

for file in 2hGlu.txt FG.txt HbA1c.txt ieu-b-110.txt ieu-b-38.txt ieu-b-39.txt IL16.txt
do
    sed -i '1s/$/\tPhenotype/' $file    # 在第一行末尾添加列名
    sed -i "2,\$s/$/\t${file%.*}/" $file    # 在第二行到最后一行末尾添加对应文件的名称
done 
```

# /share2/pub/zhenggw/zhenggw/MVMR
library(TwoSampleMR)
library(tidyverse)

GMpos <- "/share2/pub/zhenggw/zhenggw/MVMR/GM_data/"
Riskpos <- "/share2/pub/zhenggw/zhenggw/MVMR/Risk_data/"
GMfile <- c("f__Bacteroidales_noname.txt","s__Collinsella_aerofaciens.txt","s__Burkholderiales_bacterium_1_1_47.txt","s__Streptococcus_salivarius.txt","g__Burkholderiales_noname.txt","g__Collinsella.txt","s__Ruminococcus_torques.txt")
Riskfile <- c("2hGlu.txt","FG.txt","HbA1c.txt","ieu-b-110.txt","ieu-b-38.txt","ieu-b-39.txt","IL16.txt")
source("/share2/pub/zhenggw/zhenggw/MVMR/mv_extract_exposures_local.R")

for(i in GMfile){
	for(j in Riskfile){
		GMdata <- paste0(GMpos,i)
		Riskdata <- paste0(Riskpos,j)
		GM <- str_replace(i,".txt","")
		Risk <- str_replace(j,".txt","")
		fn <- paste(GM,Risk,"exposureforMVMR.txt",sep="-")

		filenames_exposure <- c(GMdata,Riskdata)
		exposure_dat <- mv_extract_exposures_local(
		  filenames_exposure,
		  sep = "\t",
		  phenotype_col = "Phenotype",
		  snp_col = "SNP",
		  beta_col = "beta",
		  se_col = "se",
		  eaf_col = "eaf",
		  effect_allele_col = "effect_allele",
		  other_allele_col = "other_allele",
		  pval_col = "pval",
		  samplesize_col = "samplesize",
		  min_pval = 1e-200,
		  log_pval = FALSE,
		  pval_threshold = 1e-05,
		  clump_r2 = 0.1,
		  clump_kb = 1000,
		  harmonise_strictness = 2
		)
	write.table(exposure_dat, file=fn, row=FALSE, col=TRUE, qu=FALSE, sep="\t")
	}
}

rm(list=ls())
expofile <- list.files(pattern="exposureforMVMR.txt$")
outlist <- c("finn-b-H7_RETINOPATHYDIAB",
			"finn-b-DM_RETINOPATHY_EXMORE",
			"finn-b-DM_RETINOPATHY")
exponame <- str_replace(expofile,"-exposureforMVMR.txt","")
resname <- paste0(exponame,"-2-",outlist,".txt")

res_list <- list() # 创建一个空列表
for(i in expofile){
  exposure_dat <- read.table(i,header=T)
  for(j in outlist){
  	exponame <- str_replace(i,"-exposureforMVMR.txt","")
	fn <- paste0("./MVMRresult/",exponame,"-2-",j,".txt")

    outcome_dat <- extract_outcome_data(exposure_dat$SNP, j)
    mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
    res <- mv_multiple(mvdat,pval_threshold = 1e-05)
    res$result$id.exposure <- exponame

    write.table(res, file=fn, sep="\t", quote=F, row.names=F)
    if(!is.null(res)){ # 检查res是否为空
      res_list <- append(res_list, res) # 将非空的res添加到列表中
    }
  }
}

# 此外，为避免在出现错误或意外情况时意外保存空结果，可以添加一个检查来判断res_list是否为空。如果不为空，则将结果保存到文件中。例如：

if(length(res_list) > 0){ # 检查res_list是否为空
  merged_res <- do.call(rbind, res_list)
  write.table(merged_res, "merged_results.txt", sep="\t", quote=F, row.names=F)
} else {
  print("No results to save") # 提示没有结果可保存
}
