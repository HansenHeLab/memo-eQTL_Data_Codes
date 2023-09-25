#############################################
## significant memo-eQTLs and characterizing
## generating figures for Fig2 and Fig2S
rm(list = ls())
setwd("./Results")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(gplots)

## readin all memo-eQTL mapping results 
res <- readRDS("./sig_memo-eQTLs.RDS")

########################################
### significant memo-eQTLs visualization
########################################
{
  ###### function ######
  plot_memo_eqtl <- function(snp_geno, cpg_meth, gene_expr, effs, pvals, hl_cut, name)
  {
    # snp_geno:  snp genotype info
    # cpg_me: CpG site methylative level 
    # gene_expr: Gene expression level       
    # res_pval : matrixeQTL p_value for Hihg_5mC a and Low_5mC group 
    # name: plot name prefix
    
    library(ggplot2)
    library(RColorBrewer)
    
    
    name_split <- strsplit(name, "_")
    name_o <- sub(":", "-", name)   ## replacing : to be able to sync to OneDrive
    
    name_e <- paste(name_split[[1]][1], name_split[[1]][3], sep = "_")    ## for eQTL
    name_me <- paste(name_split[[1]][1], name_split[[1]][2], sep = "_")   ## for meQTL
    
    
    ## split to 5mC high and low group 
    ## rmove sample equal to median cpg_meth
    #idx_mh <- cpg_meth > median(cpg_meth)
    
    ##optimal cut off
    idx_mh <- cpg_meth >= hl_cut
    idx_mh[is.na(idx_mh)] <- FALSE 
    
    idx_ml <- cpg_meth < median(cpg_meth)
    idx_ml[is.na(idx_ml)] <- FALSE 
    
    ## expr range to make sure using equal y axis
    ymax <- max(log2(gene_expr + 1)) 
    ymin <- min(log2(gene_expr + 1)) 
    
    y_range <- ymax - ymin
    ymax <- ymax + y_range * 0.05
    ycnt <- ymin - y_range * 0.05       # geno_cnt y position
    ymin <- ymin - y_range * 0.1
    
    #####################
    ## for High 5mC group 
    dat <- data.frame(geno = snp_geno[idx_mh], trait = log2(gene_expr[idx_mh] + 1))
    
    ## number of sample for each genotyp 
    geno_cnt <- table(dat$geno)
    
    gh <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
    gh <- gh + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
    gh <- gh + labs(title = name, y = "log2(FPKM +1)", x = "SNP Genotype",
                    subtitle = paste0("relH: mCpG >= ", hl_cut, "; EZ = ", effs[1], "; Pval = ", pvals[1]))
    gh <- gh + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
    gh <- gh + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
    gh <- gh + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
    gh <- gh + theme_bw() +  theme(legend.position = "none")
    
    #####################
    ## for Low 5mC group 
    dat <- data.frame(geno = snp_geno[idx_ml], trait = log2(gene_expr[idx_ml] + 1))
    
    ## number of sample for each genotyp 
    geno_cnt <- table(dat$geno)
    
    gl <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
    gl <- gl + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
    gl <- gl + labs(title =  name, y = "log2(FPKM +1)", x = "SNP Genotype",
                    subtitle = paste0("relL: mCpG < ", hl_cut, "; EZ = ", effs[2], "; Pval = ", pvals[2]))
    gl <- gl + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
    gl <- gl + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
    gl <- gl + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
    gl <- gl + theme_bw() +  theme(legend.position = "none")
    
    
    #####################
    ## for eQTL
    dat <- data.frame(geno = snp_geno, trait = log2(t(gene_expr) + 1))
    colnames(dat) <- c("geno", "trait")
    
    ## number of sample for each genotyp 
    geno_cnt <- table(dat$geno)
    
    ge <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
    ge <- ge + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
    ge <- ge + labs(title = name_e, y = "log2(FPKM +1)", x = "SNP Genotype",
                    subtitle = paste0("eQTL: Effect size = ", effs[3], "; P_val = ", pvals[3]))
    ge <- ge + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
    ge <- ge + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
    ge <- ge + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
    ge <- ge + theme_bw() +  theme(legend.position = "none")
    
    ####################
    ## for meQTL
    dat <- data.frame(geno = snp_geno, trait = t(cpg_meth))
    colnames(dat) <- c("geno", "trait")
    
    
    ## number of sample for each genotyp 
    geno_cnt <- table(dat$geno)
    
    ## beta range to make sure using equal y axis
    ymax <- max(cpg_meth, na.rm = T) 
    ymin <- min(cpg_meth, na.rm = T) 
    
    y_range <- ymax - ymin
    ymax <- ymax + y_range * 0.05
    ycnt <- ymin - y_range * 0.05       # geno_cnt y position
    ymin <- ymin - y_range * 0.1
    
    
    gme <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
    gme <- gme + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
    gme <- gme + labs(title = name_me, y = "CpG Methylation Beta value", x = "SNP Genotype",
                      subtitle = paste0("meQTL: Effect size = ", effs[4], "; P_val = ", pvals[4]))
    gme <- gme + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
    gme <- gme + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
    gme <- gme + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
    gme <- gme + theme_bw() +  theme(legend.position = "none")
    
    ## put them together 
    g_pool <- ggarrange(gh, gl, ge, gme, ncol = 2, nrow = 2)
    ggsave(paste0(name_o, "_memo-eQTLs_eQTL_meQTL.pdf"), width = 8, height = 8)
    
    ## 
    g_pool1 <- ggarrange(gh, gl, ncol = 2, nrow = 1)
    ggsave(paste0(name_o, "_memo-eQTLs.pdf"), width = 8, height = 4)
    
    g_pool2 <- ggarrange(ge, gme, ncol = 2, nrow = 1)
    ggsave(paste0(name_o, "_eQTL_meQTL.pdf"), width = 8, height = 4)
    
  }
  
  file_edata <- "./CPGEA_edata_normal_filtered.txt"
  file_gdata <- "./CPGEA_gdata_rsID_filtered.txt"
  file_mdata <- "./CPGEA_mdata_normal_filtered.txt"
  file_snp_loci <- "./CPGEA_gdata_rsID_filtered_LD_0.8_pruned_in_ATAC.bed"
  
  edata <- read.table(file_edata, header = T)
  gdata <- read.table(file_gdata, header = T)
  mdata <- read.table(file_mdata, header = T)
  snp_bed <- read.table(file_snp_loci, as.is = T)
  
  ##################################################################
  ## Draw the memo-eQTLs shown in Fig2 and and Fig 3S
  ##################################################################
  snp_cpg_gene  <- c("rs28452766_chr8:98471622_OSR2",
                     "rs1932649_chr9:84238845_IDNK",
                     "rs1507970_chr4:54683814_SRD5A3") 
       
  L <- length(snp_cpg_gene)
  
  for(i in 1:L)
  {
    ## res from matrixeQTL
    name <- snp_cpg_gene[i]    ## for me-eqtl
    
    idx_s <- match(snp_cpg_gene[i], paste(res$svar,  res$mvar,  res$evar, sep= "_"))
    
    effs <- c(round(res$clus2.coef.svar.L[idx_s], 2), 
              round(res$clus1.coef.svar.L[idx_s], 2),
              round(res$m1_svar_beta[idx_s], 2),
              round(res$m4svarbeta[idx_s], 2))
    
    pvals <- c(formatC(res$clus2.pval[idx_s], format = "e", digits = 2), 
               formatC(res$clus1.pval[idx_s], format = "e", digits = 2),
               formatC(res$m1_model_pval[idx_s], format = "e", digits = 2),
               formatC(res$m4svarpval[idx_s], format = "e", digits = 2))
    
    hl_cut <- as.numeric(strsplit(res$clus2label[idx_s], "_")[[1]][2])
    
    ## IDs
    snp_id <- strsplit(snp_cpg_gene[i], "_")[[1]][1]
    cpg_id <- strsplit(snp_cpg_gene[i], "_")[[1]][2]
    gene_id <- strsplit(snp_cpg_gene[i], "_")[[1]][3]
    
    ## to genotyp 
    idx_snp <- match(snp_id, rownames(gdata))
    idx_snp_bed <- match(snp_id, snp_bed$V4)
    homo_ref <- paste0(snp_bed$V7[idx_snp_bed], snp_bed$V7[idx_snp_bed]) 
    homo_alt <- paste0(snp_bed$V8[idx_snp_bed], snp_bed$V8[idx_snp_bed])
    hete  <- paste0(snp_bed$V7[idx_snp_bed], snp_bed$V8[idx_snp_bed])
    
    #snp_geno <- factor(gdata[idx_snp, ], labels = c(homo_ref, hete, homo_alt))'
    ## could be only two levels 
    snp_geno <- factor(gdata[idx_snp, ])
    levels(snp_geno)[levels(snp_geno) == "0"] <- homo_ref
    levels(snp_geno)[levels(snp_geno) == "1"] <- hete
    levels(snp_geno)[levels(snp_geno) == "2"] <- homo_alt
    
    idx_cpg <- match(cpg_id, rownames(mdata))
    cpg_meth <- as.matrix(mdata[idx_cpg, ])
    
    idx_gene <- match(gene_id, rownames(edata))
    gene_expr <- as.matrix(edata[idx_gene, ])
    
    plot_memo_eqtl(snp_geno, cpg_meth, gene_expr, effs, pvals, hl_cut, name)
  }
  
}


