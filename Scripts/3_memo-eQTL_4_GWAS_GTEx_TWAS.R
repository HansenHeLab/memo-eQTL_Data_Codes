########################################################
## checking memo-eQTLs in GWAS and GTEX prostate results
## and TWAS analysis 
## generating figures for Fig3A

rm(list = ls())
setwd("./Results")

library(dplyr)
library(ggpubr)

## readin all memo-eQTL mapping results 
res <- readRDS("./sig_memo-eQTLs.RDS")

############################################
## count of eSNP, eCpG and eGene
###########################################
{
  cnt_1 <- table(res$evar)
  idx_1 <- cnt_1 == 1
  tt <- paste0( sum(!idx_1), " of ", 
               length(cnt_1), " unique eGene reported >= 2")
  pdf("eGene_greater_than2_histplot.pdf", width = 5, height = 3)
  hist(cnt_1[!idx_1], main = tt, cex.main = 0.75, xlab = "# ccurence" )
  dev.off()
  
  pdf("eGene_greater_than2_histplot_pub.pdf", width = 4, height = 3)
  hist(cnt_1[!idx_1], main = "", cex.main = 0.75, xlab = "Occurence of eGene", col = "skyblue" )
  dev.off()
  
  max(cnt_1)
  
  cnt_2 <- table(res$svar)
  idx_2 <- cnt_2 == 1
  tt <- paste0( sum(!idx_2), " of ", 
                length(cnt_2), " unique eSNP reported >= 2")
  pdf("eSNP_greater_than2_histplot.pdf", width = 5, height = 3)
  hist(cnt_2[!idx_2], main = tt, cex.main = 0.75, xlab = "# ccurence" )
  dev.off()
  
  pdf("eSNP_greater_than2_histplot_pub.pdf", width = 4, height = 3)
  hist(cnt_2[!idx_2], main = "", cex.main = 0.75, xlab = "Occurence of eSNP", col = "skyblue" )
  dev.off()
  
  max(cnt_2)
  
  cnt_3 <- table(res$mvar)
  idx_3 <- cnt_3 == 1
  tt <- paste0( sum(!idx_3), " of ", 
                length(cnt_3), " unique eCpG reported >= 2")
  pdf("eCpG_greater_than2_histplot.pdf", width = 5, height = 3)
  hist(cnt_3[!idx_3], main = tt, cex.main = 0.75, xlab = "# ccurence" )
  dev.off()
  
  pdf("eCpG_greater_than2_histplot_pub.pdf", width = 4, height = 3)
  hist(cnt_3[!idx_3], main = "", cex.main = 0.75, xlab = "Occurence of eCpG", col = "skyblue" )
  dev.off()
  
  max(cnt_3)
 
}


##################################
## overlapping with all GWAS rSNPs
##################################
{
  gwas <- read.csv("../Public_Data/GWAS/GWAS_all_associations_v1.0_light_rsID_TRAIT_SNP_pairs.csv", header = F)
  colnames(gwas) <- c("trait", "snp")
  snp_trait <- paste0(gwas$snp, gwas$trait)
  length(unique(snp_trait))
  
  dim(gwas)
  length(unique(gwas$snp)) 
  length(unique(gwas$trait)) 
  
  idx <- match(gwas$snp, res$svar, nomatch = 0)
  gwas_me <- gwas[idx, ]
  write.csv(gwas_me, "memo-eQTLs_ovlapped_with_GWAS.csv")
  
  dim(gwas_me)
  length(unique(gwas_me$snp)) 
  length(unique(gwas_me$trait)) 
  snp_m_gwas <- unique(gwas_me$snp)     ## matched SNPs
  
  
  ## LD SNP in ASN 
  gwas_ld <- read.table("../Public_Data/GWAS/merged_GWAS_LD_SNPs_EAS_hg38_rsID.txt")
  length(unique(gwas_ld$V1))
  
  idx_ld <- match(gwas_ld$V1, res$svar, nomatch = 0)
  gwas_me_ld <- gwas_ld[idx_ld, ]
  
  length(gwas_me_ld)
  length(unique(gwas_me_ld)) 
  
  ### PCa risk SNPs
  gwas_pca <- read.table("../Public_Data/GWAS/PCa_269_rSNPs_hg38.bed", as.is = T)
  pca_rsnp <- gwas_pca$V4
  
  idx <- match(pca_rsnp, res$svar, nomatch = 0)
  pca_rsnp[idx]
}


##############################################
## overlapping with original GTEX PCa sig eQTL
## Prostate only
##############################################
{
  
  gtex_sig_pair <- read.table("../Public_Data/GTEx/Prostate.v8.signif_variant_gene_pairs_ID_only.txt", header = T)
  
  egene <- read.table("../Public_Data/GTEx/Prostate.v8.egenes.txt", header = T, sep = "\t")
  
  idx_g <- match(gtex_sig_pair$gene_id, egene$gene_id)
  gene_name <- egene$gene_name[idx_g]
  gtex_sig_pair <- cbind(gtex_sig_pair, gene_name)
  
  gtex_snp_gene <- paste(gtex_sig_pair$variant_id, gtex_sig_pair$gene_name, sep = "_")
  
  ## id to match GTEX SNP ID
  res_snp <- paste(res$chr, res$snp_pos, res$snp_ref, res$snp_alt, "b38", sep = "_")
  res_snp_gene <- paste(res$chr, res$snp_pos, res$snp_ref, res$snp_alt, "b38", res$evar, sep = "_")
  
  ## matched SNP 
  idx <- match(res_snp, gtex_sig_pair$variant_id)
  sum(is.na(idx))
  ## unique SNP 
  length(unique(res_snp[is.na(idx)]))
  
  ## match SNP-gene pair
  idx <- match(res_snp_gene, gtex_snp_gene)
  sum(is.na(idx))
  
  ###
  
  gtex_sig_pair <- read.table("../Public_Data/GTEx/Prostate.v8.signif_variant_gene_pairs_ID_only.txt", header = T)
  
  egene <- read.table("../Public_Data/GTEx/Prostate.v8.egenes.txt", header = T, sep = "\t")
  
  idx_g <- match(gtex_sig_pair$gene_id, egene$gene_id)
  gene_name <- egene$gene_name[idx_g]
  gtex_sig_pair <- cbind(gtex_sig_pair, gene_name)
  
  gtex_snp_gene <- paste(gtex_sig_pair$variant_id, gtex_sig_pair$gene_name, sep = "_")
  
  ## id to match GTEX SNP ID
  res_snp <- paste(res$chr, res$snp_pos, res$snp_ref, res$snp_alt, "b38", sep = "_")
  res_snp_gene <- paste(res$chr, res$snp_pos, res$snp_ref, res$snp_alt, "b38", res$evar, sep = "_")
  
  ## non-matched SNP 
  idx <- match(res_snp, gtex_sig_pair$variant_id)
  sum(is.na(idx))
  ## unique SNP 
  length(unique(res_snp[is.na(idx)]))
  
  ## matched eGene
  idx <- match(res$evar, gtex_sig_pair$gene_name)
  sum(!is.na(idx))
  length(unique(res$evar[!is.na(idx)]))
  
  egene_m <- unique(res$evar[!is.na(idx)])
  write.csv(egene_m, "matched_eGenes_for_memo-eQTL_and_GTEx_eQTL.csv")
  # corresponding SNPs
  length(unique(res$svar[!is.na(idx)]))
  
  
  ## non-match SNP-gene pair
  idx <- match(res_snp_gene, gtex_snp_gene)
  sum(is.na(idx))
  
  ## eQTL or meQTL with p
  sum(res$m1_model_pval < 0.05)
  sum(res$m4svarpval < 0.05)
  
}

##############################################
## TWAS analysis to infer "risk genes" by 
## combining eQTL and GWAS results
##############################################
{
  ## matched_eGenes_for_memo-eQTL_and_GTEx_eQTL
  egene_m 
  memo_egene <- unique(res$evar)
  
  ## 
  s <- read.table("../Public_Data/TWAS/GCST90102439_gwas_4twas_all_res.txt", header = T)
  idx_g <- match(s$ID, egene_m)
  idx_sig <- s$TWAS.P < 0.05
  sum(idx_sig, na.rm = T)
  s[!is.na(idx_g) & idx_sig , ] 
  
  idx_memo_egene <- match(s$ID, memo_egene)
  sum(!is.na(idx_memo_egene) & idx_sig, na.rm = T)
  
  ## 
  s <- read.table("../Public_Data/TWAS/GCST90011808_gwas_ldsc_all_res.txt", header = T)
  idx_g <- match(s$ID, egene_m)
  idx_sig <- s$TWAS.P < 0.05
  sum(idx_sig, na.rm = T)
  s[!is.na(idx_g) & idx_sig , ] 
  
  idx_memo_egene <- match(s$ID, memo_egene)
  sum(!is.na(idx_memo_egene) & idx_sig, na.rm = T)
  

  
  
}
