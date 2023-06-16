mv_extract_exposures_local <- function(filenames_exposure, sep = " ", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", units_col = "units", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize", gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, pval_threshold=5e-8, clump_r2=0.001, clump_kb=10000, harmonise_strictness=2)
{
	message("WARNING: Experimental function")
	l_full <- list()
	l_inst <- list()
	for(i in 1:length(filenames_exposure))
	{
		l_full[[i]] <- read_outcome_data(filenames_exposure[i], 
			sep = sep,
			phenotype_col = phenotype_col,
			snp_col = snp_col,
			beta_col = beta_col,
			se_col = se_col,
			eaf_col = eaf_col,
			effect_allele_col = effect_allele_col,
			other_allele_col = other_allele_col,
			pval_col = pval_col,
			samplesize_col = samplesize_col,
		)
		l_inst[[i]] <- subset(l_full[[i]], pval.outcome < pval_threshold)
		l_inst[[i]] <- convert_outcome_to_exposure(l_inst[[i]])
		l_inst[[i]] <- subset(l_inst[[i]], pval.exposure < pval_threshold)
		# l_inst[[i]] <- clump_data(l_inst[[i]], clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
		temp_z <- ld_clump(dplyr::tibble(rsid=l_inst[[i]]$SNP, pval=l_inst[[i]]$pval.exposure),clump_p=pval_threshold,
			    plink_bin = "/share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps/plink",
			    bfile = "/share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps/EUR_ref/EUR",
			    clump_kb = clump_kb,clump_r2 = clump_r2)
	    l_inst[[i]] <- l_inst[[i]][which(l_inst[[i]]$SNP %in% temp_z$rsid),]
	}

	exposure_dat <- dplyr::bind_rows(l_inst)
	id_exposure <- unique(exposure_dat$id.exposure)
	temp <- exposure_dat
	temp$id.exposure <- 1
	temp <- temp[order(temp$pval.exposure, decreasing=FALSE), ]
	temp <- subset(temp, !duplicated(SNP))
	temp_z <- ld_clump(dplyr::tibble(rsid=temp$SNP, pval=temp$pval.exposure),clump_p=pval_threshold,
			    plink_bin = "/share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps/plink",
			    bfile = "/share2/pub/zhenggw/zhenggw/GM_GWAS_LD_clumped_snps/EUR_ref/EUR",
			    clump_kb = clump_kb,clump_r2 = clump_r2)
    temp <- temp[which(temp$SNP %in% temp_z$rsid),]
	# temp <- clump_data(temp, clump_p1=pval_threshold, clump_r2=clump_r2, clump_kb=clump_kb)
	exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)

	d1 <- lapply(l_full, function(x) {
		subset(x, SNP %in% exposure_dat$SNP)
		}) %>% dplyr::bind_rows()

	stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
	d1 <- subset(d1, mr_keep.outcome)
	d2 <- subset(d1, id.outcome != id_exposure[1])
	d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))

	# Harmonise against the first id
	d <- harmonise_data(d1, d2, action=harmonise_strictness)

	# Only keep SNPs that are present in all
	tab <- table(d$SNP)
	keepsnps <- names(tab)[tab == length(id_exposure)-1]
	d <- subset(d, SNP %in% keepsnps)
	
	# Reshape exposures
	dh1 <- subset(d, id.outcome == id.outcome[1], select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure))
	dh2 <- subset(d, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))
	names(dh2) <- gsub("outcome", "exposure", names(dh2))
	dh <- rbind(dh1, dh2)
	return(dh)
}