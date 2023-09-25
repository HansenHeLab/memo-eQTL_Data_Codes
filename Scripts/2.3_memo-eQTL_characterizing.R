###################################################
## significant memo-eQTLs and characterizing
## generating figures for Fig2 and Fig2S and Fig3BC
rm(list = ls())
setwd("./Results")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(gplots)
library(clusterProfiler)

## readin all memo-eQTL mapping results 
cmp <- readRDS("./all_tested_snp_cpg_gene_cmb_4_memo-eqtl.RDS")

##########################
### significant memo-eQTLs
##########################
{
  ## m3 needs to be significant
  idx_sig <- cmp$m3_model_pval < 0.05
  cmp <- cmp[idx_sig, ]
  
  ## require m13 > 0.05 and cmp$m2m3pval
  idx_1 <- cmp$m1m3pval >= 0.05 & cmp$m2m3pval  < 0.05
  idx_2 <- cmp$m1m3pval < 0.05 & cmp$m2m3pval  < 0.05 
  idx_3 <- cmp$m1m3pval < 0.05 & cmp$m2m3pval  >= 0.05
  idx_4 <- cmp$m1m3pval >= 0.05 & cmp$m2m3pval  >= 0.05
  
  col_en <- rep("gray80", nrow(cmp))
  col_en[idx_1] <- "gray40"
  col_en[idx_2] <- "red3"
  col_en[idx_3] <- "gray60"
  col_en[idx_4] <- "gray80"
  
  table(col_en)
  
  png(file = "M13_vs_M23_Pval_for_memo-eQTLs.png", heigh = 4.5, width = 4.5, res = 600, units = "in")
  par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
  plot(-log10(cmp$m1m3pval), -log10(cmp$m2m3pval), col = col_en, pch = 16, 
       xlab = "-log10(Pval of M3 vs M1)", ylab = "-log10(Pval of M3 vs)")
  abline(v = -log10(0.05), lty = 2 , col = "gray40") 
  abline(h = -log10(0.05), lty = 2 , col = "gray40")  
  dev.off()
  
  ## M3 R2
  col_enf <- factor(col_en, levels = c("gray80", "gray60", "gray40", "red3"))       ## convert to factors
  Group <- col_enf
  levels(Group) <- c("G4", "G3", "G2", "G1")
  
  dat <- data.frame(cmp$m3_model_rsq, -log10(cmp$m3_model_pval), 
                    cmp$m3_svar_beta, -log10(cmp$m3_svar_pval), cmp$svar_percent, 
                    cmp$m3_mvar_beta, -log10(cmp$m3_mvar_pval), cmp$mvar_percent, 
                    cmp$m3_svar_mvar_beta, -log10(cmp$m3_svar_mvar_pval), cmp$svar.mvar_percent,
                    col_enf, Group)
  ylab <- c("Squared R of M3", "-log10(Pval of M3)", 
            "Beta of SNP in M3", "-log10(Pval of SNP in M3)", "Relative Var explained by SNP",
            "Beta of meCpG in M3", "-log10(Pval of meCpG in M3)", "Relative Var explained by meCpG",
            "Beta of SNP ✕ meCpG in M3", "-log10(Pval of SNP x meCpG in M3)", "Relative Var explained by SNP x meCpG",
            "Color", "Group")
  
  my_cmp <- list( c("G4", "G3"), c("G4", "G2"), c("G4", "G1"), 
                  c("G3", "G2"), c("G3", "G1"), c("G2", "G1"))
  
  L <- ncol(dat)
  for(i in 1: (L-2))
  {
    g <- ggviolin(dat, x = "Group", y = colnames(dat)[i], fill = "Group", ylab = ylab[i],
                  palette = c("gray80", "gray60", "gray40", "red3"), add = "boxplot", add.params = list(fill = "white")) 
    g <- g + stat_compare_means(comparisons = my_cmp, label = "p.signif") + theme_bw() + theme(legend.position = "top")
    ggsave(paste0("CMP_", colnames(dat)[i], "_in_4groups.pdf"), width = 5, height = 4.5)
  }
  
  
  ########################
  ## significant memo-eQTL
  {
    idx_memo <- col_en == "red3"
    cmp_memo <- cmp[idx_memo, ]
    dim(cmp_memo)
    
    saveRDS(cmp_memo, file = "sig_memo-eQTLs.RDS")
    write.csv(cmp_memo, "sig_memo-eQTLs.csv")
  }
  
  
}


##############################################
### memo-eGenes functional enrichment analysis 
#############################################
{
  snp_s <-  unique(cmp_memo$svar)
  cpg_s <-  unique(cmp_memo$mvar)
  gene_s <- unique(cmp_memo$evar) 
  
  length(snp_s)
  length(cpg_s)
  length(gene_s)
  
  {
    
    ###################################################
    ## DEGs KEGG and GO functional enrichment analysis 
    POS_gmt <- read.gmt("../Public_Data/msigdb_v7.1_GMTs/c1.all.v7.1.symbols.gmt")
    KEGG_gmt <- read.gmt("../Public_Data/msigdb_v7.1_GMTs/c2.cp.kegg.v7.1.symbols.gmt")
    BP_gmt <- read.gmt("../Public_Data/msigdb_v7.1_GMTs/c5.bp.v7.1.symbols.gmt")
    CC_gmt <- read.gmt("../Public_Data/msigdb_v7.1_GMTs/c5.cc.v7.1.symbols.gmt")
    MF_gmt <- read.gmt("../Public_Data/msigdb_v7.1_GMTs/c5.mf.v7.1.symbols.gmt")
    
    ### all eGenes
    prefix_n = "All_memo-eGene"
    {
      ##########################
      ## for memo-egenes
      pos_enrich <-  enricher(gene_s, TERM2GENE = POS_gmt) 
      if(length(pos_enrich)){
        g <- dotplot(pos_enrich) + ggtitle(paste0(prefix_n, "_POS"))
        ggsave(paste0(prefix_n, "_POS.pdf"), width = 5, height = 4)
      }
      
      kegg_enrich <-  enricher(gene_s, TERM2GENE = KEGG_gmt) 
      if(length(kegg_enrich$geneID)){
        g <- dotplot(kegg_enrich) + ggtitle(paste0(prefix_n, "_KEGG"))
        ggsave(paste0(prefix_n, "_KEGG.pdf"), width = 9, height = 4)
      }
      
      bp_enrich <-  enricher(gene_s, TERM2GENE = BP_gmt) 
      if(length(bp_enrich$geneID)){
        g <- dotplot(bp_enrich) + ggtitle(paste0(prefix_n, "_GO_BP"))
        ggsave(paste0(prefix_n, "_GO_BP.pdf"), width = 11, height = 4)
      }
      
      cc_enrich <-  enricher(gene_s, TERM2GENE = CC_gmt) 
      if(length(cc_enrich$geneID)){
        g <- dotplot(cc_enrich) + ggtitle(paste0(prefix_n, "_GO_CC"))
        ggsave(paste0(prefix_n, "_GO_CC.pdf"), width = 9, height = 4)
      }
      
      mf_enrich <-  enricher(gene_s, TERM2GENE = MF_gmt) 
      if(length(mf_enrich$geneID)){
        g <- dotplot(mf_enrich) + ggtitle(paste0(prefix_n, "_GO_MF"))
        ggsave(paste0(prefix_n, "_GO_MF.pdf"), width = 9, height = 4)
      }
      
    }
    
    
    ### for the eGenes enriched in top2 chromatin loci
    {
      ## for the gene enriched in chr6p21
      gene_pos <- unlist(strsplit(pos_enrich$geneID[1], "\\/"))      ## chr6p21
      gene_pos_kegg <-  enricher(gene_pos, TERM2GENE = KEGG_gmt) 
      
      if(length(gene_pos_kegg$geneID)){
        g <- dotplot(gene_pos_kegg) + ggtitle("chr6p21_eGenes_KEGG_enrichment")
        ggsave("chr6p21_me_eGenes_KEGG_enrichment.pdf", width = 7, height = 4)
      }
      
      gene_pos_bp <-  enricher(gene_pos, TERM2GENE = BP_gmt) 
      if(length(gene_pos_bp$geneID)){
        g <- dotplot(gene_pos_bp) + ggtitle("chr6p21_eGenes_BP_enrichment")
        ggsave("chr6p21_me_eGenes_BP_enrichment.pdf", width = 7, height = 4)
      }
      
      ## for the gene enriched in chr17q25  !! no enrichment
      gene_pos <- unlist(strsplit(pos_enrich$geneID[2], "\\/"))      ## chr6p21
      gene_pos_kegg <-  enricher(gene_pos, TERM2GENE = KEGG_gmt) 
      
      if(length(gene_pos_kegg$geneID)){
        g <- dotplot(gene_pos_kegg) + ggtitle("chr17q25_eGenes_KEGG_enrichment")
        ggsave("chr17q25_me_eGenes_KEGG_enrichment.pdf", width = 7, height = 4)
      }
      
      gene_pos_bp <-  enricher(gene_pos, TERM2GENE = BP_gmt) 
      if(length(gene_pos_bp$geneID)){
        g <- dotplot(gene_pos_bp) + ggtitle("cchr17q25_eGenes_BP_enrichment")
        ggsave("chr17q25_me_eGenes_BP_enrichment.pdf", width = 7, height = 4)
      }
    }
    
    ### for rest genes other than those located in chr6p21
    {
      gene_rm <- unlist(strsplit(pos_enrich$geneID[1], "\\/")) 
      idx_rm <- match(gene_s, gene_rm)
      gene_ss <- gene_s[is.na(idx_rm)]
      
      prefix_n = "All_memo-eGene_except_chr6P21"
      {
        ##########################
        ## for memo-egenes
        pos_enrich <-  enricher(gene_ss, TERM2GENE = POS_gmt) 
        if(length(pos_enrich)){
          g <- dotplot(pos_enrich) + ggtitle(paste0(prefix_n, "_POS"))
          ggsave(paste0(prefix_n, "_POS.pdf"), width = 5, height = 4)
        }
        
        
        kegg_enrich <-  enricher(gene_ss, TERM2GENE = KEGG_gmt) 
        if(length(kegg_enrich$geneID)){
          g <- dotplot(kegg_enrich) + ggtitle(paste0(prefix_n, "_KEGG"))
          ggsave(paste0(prefix_n, "_KEGG.pdf"), width = 9, height = 4)
        }
        
        bp_enrich <-  enricher(gene_ss, TERM2GENE = BP_gmt) 
        if(length(bp_enrich$geneID)){
          g <- dotplot(bp_enrich) + ggtitle(paste0(prefix_n, "_GO_BP"))
          ggsave(paste0(prefix_n, "_GO_BP.pdf"), width = 11, height = 4)
        }
        
        cc_enrich <-  enricher(gene_ss, TERM2GENE = CC_gmt) 
        if(length(cc_enrich$geneID)){
          g <- dotplot(cc_enrich) + ggtitle(paste0(prefix_n, "_GO_CC"))
          ggsave(paste0(prefix_n, "_GO_CC.pdf"), width = 9, height = 4)
        }
        
        mf_enrich <-  enricher(gene_ss, TERM2GENE = MF_gmt) 
        if(length(mf_enrich$geneID)){
          g <- dotplot(mf_enrich) + ggtitle(paste0(prefix_n, "_GO_MF"))
          ggsave(paste0(prefix_n, "_GO_MF.pdf"), width = 9, height = 4)
        }
        
        
      }
      
    }
    
  }
  
}

