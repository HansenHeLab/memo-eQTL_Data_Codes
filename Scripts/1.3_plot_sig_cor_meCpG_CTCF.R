rm(list = ls())
setwd("./Results")

library(dplyr)
library(ggrepel)

###############
# plot function 
###############
corr_scatter_ci <- function(x, y, cor_m,x_lab, y_lab, file_name)
{
  ## ggpubr::stat_cor can be an alternative way for adding r and p
  cor <- cor.test(x , y, method = cor_m)
  pval <- cor$p.value
  if(pval > 2.2e-16){
    pval <- formatC(pval, format = "e", digits = 2)
    sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")
  } else {
    sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "")
  }
  
  dat <- data.frame(x, y)
  g <- ggplot(dat, aes(x = x, y = y)) + geom_point(color='#2980B9', size = 2)
  g <- g + geom_smooth(method=lm, formula = y~x, color='#2C3E50') + theme_classic()
  g <- g + labs(x = x_lab, y = y_lab, title = sub_t) + theme(plot.title = element_text(size = 9))
  g <- g + geom_rug(color = "gray40", lwd = 0.35)
  ggsave(file_name, width = 3, height = 3)
}

###############
## readin files
###############
{
ctcf <- read.table(gzfile("../Public_Data/ENCODE/ENCODE_ctcf_matrix.txt.gz"), header = T, row.names = 1)
meth <- read.table(gzfile("../Public_Data/ENCODE/ENCODE_methylation_matrix.txt.gz"),  header = T)
meth_cnt <- read.table(gzfile("../Public_Data/ENCODE/ENCODE_methylation_totalRead_matrix.txt.gz"),  header = T)

## read in sig cor meCpG-CTCF 
sig_cor <- read.table("./cor_res_o20_sig_filtered_mostSig_cpg_perCTCF.txt", header = T)
}

##############################################
## for significantly negative correlated pairs
##############################################
{
sig_neg <- filter(sig_cor, cor_cpg_ctcf_consis == -1 ) %>% 
           arrange(qval_global, desc(pcc_2))

## PCC_2 distribution 
pdf("sig_neg_pcc_2_distribution.pdf", width = 4, height = 3)
hist(sig_neg$pcc_2, breaks = 10, main = "cor(CTCF, meth)_Neg", xlab = "PCC")
dev.off()

## top 5 significantly pairs
tmp1 <- sig_neg[1:5, ]

## top 5 around the mean PCC_2 and cpg locate near to ctcf mid
idx_t <- sig_neg$pcc_2 > -0.7 & abs(sig_neg$cpg2ctcf_mid) < 10
sig_neg_t <- sig_neg[idx_t, ]
tmp2 <- sig_neg_t[1:5, ]

sig_plot <- rbind(tmp1, tmp2)

for(i in 1:nrow(sig_plot))
{
idx_1 <- which(rownames(ctcf) == sig_plot$ctcf_id[i])
idx_2 <- which(rownames(meth) == sig_plot$cpg[i])
idx_3 <- which(rownames(meth_cnt) == sig_plot$cpg[i])
idx_s <- meth_cnt[idx_3, ] > 20

x = as.numeric(ctcf[idx_1, idx_s ])
y = as.numeric(meth[idx_2, idx_s])
x_lab = paste0(sig_plot$ctcf_id[i], "_", sig_plot$chr[i], ":", sig_plot$ctcf_star[i], "-",  sig_plot$ctcf_end[i])
y_lab = paste0("CpG_", sig_plot$cpg[i], "_methylation_level")
sub_t <- paste0("cpg2ctcf_mid:", sig_plot$cpg2ctcf_mid[i]) 
sample = colnames(ctcf)[idx_s]

dat <- data.frame(x, y, sample)
g <- ggplot(dat, aes(x = x, y = y, label = sample)) + geom_point(color = "#2980B9", size = 2)
g <- g + geom_text_repel(size = 2)
g <- g + geom_smooth(method = lm, formula = y ~ x, color = "#2C3E50") + theme_classic()
g <- g + labs(x = x_lab, y = y_lab, title = sub_t) 
g <- g + geom_rug(color = "gray40", lwd = 0.35)
ggsave(paste0("Neg_selected_eg_",i, "labeled.pdf"), width = 5, height = 5)
}

}


##############################################
## for significantly positive correlated pairs
##############################################
{
  sig_pos <- filter(sig_cor, cor_cpg_ctcf_consis == 1 ) %>% 
    arrange(qval_global, desc(pcc_2))
  
  ## PCC_2 distribution 
  pdf("POS_sig_pcc_2_distribution.pdf", width = 4, height = 3)
  hist(sig_pos$pcc_2, breaks = 10, main = "cor(CTCF, meth)_Pos", xlab = "PCC")
  dev.off()
  
  ## top 5 significant 
  tmp1 <- sig_pos[1:5, ]
  
  ## top 5 around the mean PCC_2 and cpg locate near to ctcf mid
  summary(sig_pos$pcc_2)
  idx_t <- sig_pos$pcc_2 < 0.77 & abs(sig_pos$cpg2ctcf_mid) < 10
  sig_pos_t <- sig_pos[idx_t, ]
  tmp2 <- sig_pos_t[1:5, ]
  
  sig_plot <- rbind(tmp1, tmp2)
  
  ## draw the top 5 of significant candidates
  for(i in 1:nrow(sig_plot))
  {
    idx_1 <- which(rownames(ctcf) == sig_plot$ctcf_id[i])
    idx_2 <- which(rownames(meth) == sig_plot$cpg[i])
    idx_3 <- which(rownames(meth_cnt) == sig_plot$cpg[i])
    idx_s <- meth_cnt[idx_3, ] > 20
    
    x = as.numeric(ctcf[idx_1, idx_s ])
    y = as.numeric(meth[idx_2, idx_s])
    x_lab = paste0(sig_plot$ctcf_id[i], "_", sig_plot$chr[i], ":", sig_plot$ctcf_star[i], "-",  sig_plot$ctcf_end[i])
    y_lab = paste0("CpG_", sig_plot$cpg[i], "_methylation_level")
    sub_t <- paste0("cpg2ctcf_mid:", sig_plot$cpg2ctcf_mid[i]) 
    sample = colnames(ctcf)[idx_s]
    
    dat <- data.frame(x, y, sample)
    g <- ggplot(dat, aes(x = x, y = y, label = sample)) + geom_point(color = "#2980B9", size = 2)
    g <- g + geom_text_repel(size = 2)
    g <- g + geom_smooth(method = lm, formula = y ~ x, color = "#2C3E50") + theme_classic()
    g <- g + labs(x = x_lab, y = y_lab, title = sub_t) 
    g <- g + geom_rug(color = "gray40", lwd = 0.35)
    ggsave(paste0("Pos_selected_eg_",i, "labeled.pdf"), width = 5, height = 5)
  }
}


