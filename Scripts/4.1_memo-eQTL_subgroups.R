######################################################################
## significant memo-eQTLs subgrouping  and linear distance comaprisons
## generating figures for Fig4 and Fig4S
rm(list = ls())
setwd("./Results")

library(ggplot2)
library(dplyr)

cmp_memo <- readRDS("./sig_memo-eQTLs.RDS")

###################################################
# add cor(cpg, ctcf)  snp and gene  information  
###################################################
{
## for cor(CpG, CTCF)
cor_cpg_ctcf <- read.table("./cor_res_o20_sig_filtered_mostSig_cpg_perCTCF.txt", header = T, as.is = T)
idx_cpg <- match(cmp_memo$mvar, cor_cpg_ctcf$cpg)

cor_cpg_ctcf_add <- cor_cpg_ctcf[idx_cpg, ] %>%  dplyr::select(ctcf_id, cor_cpg_ctcf = scc_2, cor_cpg_ctcf_fdr = qval_global, 
                                                               chr, cpg_pos, ctcf_mid)

file_tss_loci <- "./CPGEA_edata_normal_filtered_TSS.bed"
file_snp_loci <- "./CPGEA_gdata_rsID_filtered_LD_0.8_pruned_in_ATAC.bed"
file_cpg_loci <- "./CPGEA_mdata_normal_filtered_CpG.bed"

## snp info 
snp_info <- read.table(file_snp_loci, header = F, as.is = T) 
idx_snp <- match(cmp_memo$svar, snp_info$V4)
snp_info_add <- snp_info[idx_snp, ] %>%  dplyr::select(snp_pos = V3, snp_ref = V7, snp_alt = V8)


## PRAD distal ATAC-seq peak info
atac_info <- read.table(file_snp_loci, as.is = T)
idx_atac <-  match(cmp_memo$svar, atac_info$V4)
atac_info_add <- atac_info[idx_atac, ] %>%    
  mutate(snp2atac_mid = V3 - (V10 + V11)/2) %>% 
  dplyr::select(atac_id = V12, atac_start = V10, atac_end = V11, snp2atac_mid) 


## gene info
gene_info <- read.table(file_tss_loci, header = F, as.is = T) 
idx_gene <- match(cmp_memo$evar, gene_info$V4)
gene_info_add <- gene_info[idx_gene, ] %>%  dplyr::select(tss_pos = V3, gene_type = V5, gene_strand = V6)

res <- cbind(cmp_memo, cor_cpg_ctcf_add, snp_info_add, atac_info_add,  gene_info_add)

idx_pos <- res$cor_cpg_ctcf > 0
idx_neg <- res$cor_cpg_ctcf < 0
sum(idx_pos)
sum(idx_neg)

res_pos <- res[idx_pos, ]
res_neg <- res[idx_neg, ]

#write.csv(res_pos, file = "memo-eqtl_only_underlie_pos_cor_cpg_ctcf.csv")
#write.csv(res_neg, file = "memo-eqtl_only_underlie_neg_cor_cpg_ctcf.csv")

#saveRDS(res_pos, file = "memo-eqtl_only_underlie_pos_cor_cpg_ctcf.rds")
#saveRDS(res_neg, file = "memo-eqtl_only_underlie_neg_cor_cpg_ctcf.rds")

}


##################################################
##  memo-eQTLs separation with neg cor(cpg, ctcf)!!
#################################################
{
  #res <- res_neg[, -c(37:39)]     ## rm duplicated columns 
  
  col_memo_lh <- rep("gray", nrow(res))
  idx_1 <- res$clus1.modelpval.value >= 0.05 & res$clus2.modelpval.value  < 0.05
  idx_2 <- res$clus1.modelpval.value < 0.05 & res$clus2.modelpval.value  < 0.05
  idx_3 <- res$clus1.modelpval.value < 0.05 & res$clus2.modelpval.value  >= 0.05
  idx_4 <- res$clus1.modelpval.value >= 0.05 & res$clus2.modelpval.value  >= 0.05
  
  col_memo_lh [idx_1] <- "coral"
  col_memo_lh [idx_2] <- "coral4"
  col_memo_lh [idx_3] <- "cadetblue"
  col_memo_lh [idx_4] <- "gray"
  
  table(col_memo_lh)
  
  ## pavlues 
  png(file = "memo-eQTLs_meLow_vs_meHigh.png", heigh = 4, width = 4, res = 600, units = "in")
  par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
  plot(-log10(res$clus1.modelpval.value), -log10(res$clus2.modelpval.value),  pch = 16, 
       xlab = "-log10(pval for relL)", ylab = "-log10(pval for relH)", col = col_memo_lh )
  #abline(a = 0 , b = 1, lty = 2 , col = "gray40")
  abline(v = -log10(0.05), lty = 2 , col = "gray40") 
  abline(h = -log10(0.05), lty = 2 , col = "gray40")  
  dev.off()
}


##############################
## me_eqtl in me_high, me_low 
#############################
{
  # display.brewer.all()

  sigGroup <- vector()
  sigGroup[idx_1] <- "sigHigh"
  sigGroup[idx_2] <- "sigBoth"
  sigGroup[idx_3] <- "sigLow"
  sigGroup[idx_4] <- "sigNone"
  sigGroup <- as.factor(sigGroup)
  sigColor <- c("coral4", "coral", "cadetblue", "gray")
  
  res <- cbind(res, sigGroup)
  
  ################################################
  ## compare the distance features between groups
  ################################################
  {
    ## snps in same ATAC peak region may yied same me_eQTLs
    res <-  mutate(res, 
                   snp2cpg = snp_pos - cpg_pos,
                   snp2gene = snp_pos - tss_pos,
                   cpg2gene = cpg_pos - tss_pos,
                   cpg2ctcf_mid = cpg_pos - ctcf_mid)           
    
    saveRDS(res, "memo-eQTLs_meta_info_1711.RDS")
    #write.csv(res, "memo-eQTLs_meta_info_1711.csv")
    
    ## cor(CpG, CTCF) in four group 
    my_cmp <- list( c("sigBoth", "sigHigh"), c("sigBoth", "sigLow"), c("sigBoth", "sigNone"), 
                    c("sigHigh", "sigLow"), c("sigHigh", "sigNone"), c("sigLow", "sigNone"))
    
    dat_cmp <- select(res, cor_cpg_ctcf, snp2cpg, snp2gene, cpg2gene, cpg2ctcf_mid, sigGroup)
    
    L <- ncol(dat_cmp)
    
    for(i in 1: (L-1))
    {
      g <- ggviolin(dat_cmp, x = "sigGroup", y = colnames(dat_cmp)[i], fill = "sigGroup",
                    palette = sigColor, add = "boxplot", add.params = list(fill = "white")) 
      g <- g + stat_compare_means(comparisons = my_cmp, label = "p.signif") 
      g <- g + theme_bw() + theme(legend.position = "none")
      ggsave(paste0("4Groups_memo-eQTLs_", colnames(dat_cmp)[i], "_cmp.pdf"), width = 5, height = 4)
    }
    
    ## high vs low only 
    idx_hl <- sigGroup == "sigHigh" | sigGroup == "sigLow"
    dat_cmp_hl <- dat_cmp[idx_hl, ]
    my_cmp_hl  <- list(c("sigHigh", "sigLow")) 
    
    L <- ncol(dat_cmp)
    
    for(i in 1: (L-1))
    {
      g <- ggviolin(dat_cmp_hl, x = "sigGroup", y = colnames(dat_cmp_hl)[i], fill = "sigGroup",
                    palette = c("coral", "cadetblue"), add = "boxplot", add.params = list(fill = "white")) 
      g <- g + stat_compare_means(comparisons = my_cmp_hl, label = "p.signif") 
      ggsave(paste0("meHigh_vs_meLow_memo-eQTLs_", colnames(dat_cmp)[i], "_cmp.pdf"), width = 3, height = 4)
    }
    
    
    ####################
    ## distance pattern
    
    #idx_h <-res$sigGroup == "High_5mC"
    res_pos_dis <-  data.frame(ungroup(res) %>%  dplyr::select(snp2cpg:cpg2ctcf_mid, sigGroup))
    N <- ncol(res_pos_dis)
    
    ## gene on minus strand 
    group_by(res, sigGroup, gene_strand) %>%  summarise(gene_cnt = n())     ## checking gene on plus and minus strand
    
    idx_minus <-res$gene_strand == "-"
    res_pos_dis_strand <-res_pos_dis
    res_pos_dis_strand[idx_minus, 1:(N-1)] <- -res_pos_dis_strand[idx_minus, 1:(N-1)]
    
    prefix <- c("Distance_from_eSNP_to_eCpG", "Distance_from_eSNP_to_eGene",
                "Distance_from_eCpG_to_eGene", "Distance_from_eCpG_to_CTCF_mid")
    
    idx_sig_high <- sigGroup == "sigHigh"
    idx_sig_low <- sigGroup == "sigLow"
    
    ## for all four groups
    for(i in 1:(N-1)){
      
      ### density of distance 
      ks <- ks.test(res_pos_dis[idx_sig_high, i],res_pos_dis[idx_sig_low, i])
      pval <- formatC(ks$p.value, format = "e", digits = 2)
      tt <- paste0("HvsL K-S test: pval = ", pval)
      g <- ggplot(res_pos_dis, aes(x =res_pos_dis[, i],  color=sigGroup)) + geom_density() 
      g <- g + scale_color_manual(values = sigColor)
      #g <- g + labs(x = paste0(prefix[i], ": < 0 leftside"), title = tt) + theme_bw() + theme(legend.position = "top") 
      g <- g + labs(x = paste0(prefix[i]), title = tt) + theme_classic() + theme(legend.position = "top") 
      ggsave(paste0(prefix[i], ".pdf"), width = 4.9, height = 4)
      
      ## relative positional sign amony snp, cgp and gene   
      ## upstream = -1, downstream = 1, same = 0 
      ### density of distance by taking gene strand into account 
      ks <- ks.test(res_pos_dis_strand[idx_sig_high, i],res_pos_dis_strand[idx_sig_low, i])
      pval <- formatC(ks$p.value, format = "e", digits = 2)
      tt <- paste0("HvsL K-S test: pval = ", pval)
      g <- ggplot(res_pos_dis_strand, aes(x =res_pos_dis_strand[, i],  color=sigGroup)) + geom_density() 
      g <- g + scale_color_manual(values = sigColor)
      g <- g + labs(x = paste0(prefix[i], ": < 0 upstream"), title = tt) + theme_bw()  + theme(legend.position = "bottom") 
      ggsave(paste0(prefix[i], "_consider_gene_strand.pdf"), width = 3.5, height = 3)
      
    }
    
    
    ## for high vs low groups 
    res_pos_dis_hl <- res_pos_dis[idx_hl, ]
    res_pos_dis_strand_hl <- res_pos_dis_strand[idx_hl, ]
    prefix <- c("meHigh_vs_meLow_Distance_from_SNP_to_CpG", "meHigh_vs_meLow_Distance_from_SNP_to_Gene",
                "meHigh_vs_meLow_Distance_from_CpG_to_Gene", "meHigh_vs_meLow_Distance_from_CpG_to_CTCF_mid")
    
    for(i in 1:(N-1)){
      
      ### density of distance 
      ks <- ks.test(res_pos_dis[idx_sig_high, i],res_pos_dis[idx_sig_low, i])
      pval <- formatC(ks$p.value, format = "e", digits = 2)
      tt <- paste0("HvsL K-S test: pval = ", pval)
      g <- ggplot(res_pos_dis_hl, aes(x =res_pos_dis_hl[, i],  color=sigGroup)) + geom_density() 
      g <- g + scale_color_manual(values = c("coral", "cadetblue"))
      g <- g + labs(x = paste0(prefix[i], ": < 0 leftside"), title = tt) + theme_bw() + theme(legend.position = "bottom") 
      ggsave(paste0(prefix[i], ".pdf"), width = 3.5, height = 3)
      
      ## relative positional sign amony snp, cgp and gene   
      ## upstream = -1, downstream = 1, same = 0 
      ### density of distance by taking gene strand into account 
      ks <- ks.test(res_pos_dis_strand[idx_sig_high, i],res_pos_dis_strand[idx_sig_low, i])
      pval <- formatC(ks$p.value, format = "e", digits = 2)
      tt <- paste0("HvsL K-S test: pval = ", pval)
      g <- ggplot(res_pos_dis_strand_hl, aes(x =res_pos_dis_strand_hl[, i],  color=sigGroup)) + geom_density() 
      g <- g + scale_color_manual(values = c("coral", "cadetblue"))
      g <- g + labs(x = paste0(prefix[i], ": < 0 upstream"), title = tt) + theme_bw()  + theme(legend.position = "bottom") 
      ggsave(paste0(prefix[i], "_consider_gene_strand.pdf"), width = 3.5, height = 3)
      
    }
    
    ## for high vs low groups for CpG_to_CTCF_mid
    prefix <- c("meHigh_vs_meLow_Distance_from_SNP_to_CpG", "meHigh_vs_meLow_Distance_from_SNP_to_Gene",
              "meHigh_vs_meLow_Distance_from_CpG_to_Gene", "Distance_from_eCpG_to_CTCF")
    for(i in 4){
      
      ### density of distance 
      ks <- ks.test(res_pos_dis[idx_sig_high, i],res_pos_dis[idx_sig_low, i])
      pval <- formatC(ks$p.value, format = "e", digits = 2)
      tt <- paste0("HvsL K-S test: pval = ", pval)
      tt
      
      g <- ggplot(res_pos_dis_hl, aes(x =res_pos_dis_hl[, i],  color=sigGroup)) + geom_density() 
      g <- g + scale_color_manual(values = c("coral", "cadetblue"))
      g <- g + labs(x = paste0(prefix[i])) + theme_classic() + theme(legend.position = "top") 
      ggsave(paste0(prefix[i], ".pdf"), width = 4.5, height = 4)
      
    }

   
}
  

}

###################################
##  beta in high and low meth group
###################################
{
  ## all memo-eQTLs
  png(file = "memo-eQTLs_meLow_vs_meHigh_beta.png", heigh = 4.5, width = 4.5, res = 600, units = "in")
  par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
  plot(res$clus1.coef.svar.L, res$clus2.coef.svar.L,  pch = 16, 
       xlab = "meLow-SNP-beta)", ylab = "meHigh-SNP-beta", col = "red3")
  #abline(a = 0 , b = 1, lty = 2 , col = "gray40")
  abline(v = 0, lty = 2 , col = "gray40") 
  abline(h = 0, lty = 2 , col = "gray40")  
  dev.off()
  
  ## four groups
  idx_sig_both <- sigGroup == "sigBoth" 
  idx_sig_none <- sigGroup == "sigNone" 
  
  idx_sig <- list(idx_sig_low, idx_sig_both, idx_sig_high, idx_sig_none)
  idx_sig_name <- unique(sigGroup)
  idx_sig_col <- c("cadetblue", "coral4", "coral", "gray")
  
  for (i in 1:4)
  {
  name = paste0("memo-eQTLs_meLow_vs_meHigh_", idx_sig_name[i], "_beta.png")
  png(file = name, heigh = 4, width = 4, res = 600, units = "in")
  par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
  plot(res$clus1.coef.svar.L[idx_sig[[i]]], res$clus2.coef.svar.L[idx_sig[[i]]],  pch = 16, 
       xlab = "Effect Size for relL", ylab = "Effect Size for relH", col = idx_sig_col[i])
  #abline(a = 0 , b = 1, lty = 2 , col = "gray40")
  abline(v = 0, lty = 2 , col = "gray40") 
  abline(h = 0, lty = 2 , col = "gray40")  
  dev.off()
  }
  
}




