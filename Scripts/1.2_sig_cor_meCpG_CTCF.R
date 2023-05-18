rm(list = ls())
setwd("./Results")

library(ggplot2)
library(dplyr)
library(scales)

###############################
## pre-filtering and annotation
###############################
{
cor_res <- read.table(gzfile("./cor_meCpG_CTCF_all.txt.gz"), as.is = T)

## pcc_1 and pval_2 are pearson correlation coefficient and p vale for all meCpG-CTCF pairs; 
## pcc_2 and pval_2 are pearson correlation coefficient and p vale for meCpG-CTCF pairs 
##   with methyalation total reads >= 20 in at least 10 samples 
colnames(cor_res) <- c("ctcf_id", "cpg", "numTested_all", "pcc_1", "pval_1", "numTested_o20", "pcc_2", "pval_2")

#########################################
## Focusing on variable meCpGs and CTCFs
{
## CpG methylation levels
cpg_val <- read.table(gzfile("../Public_Data/ENCODE/ENCODE_methylation_matrix.txt.gz"), header = T, row.names = 1)
cpg_val_iqr<- apply(cpg_val, 1,  IQR, na.rm = T)

pdf(file = "meCpGs_IQR.pdf", width = 5, height = 3)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(cpg_val_iqr, main = "", xlab = "IQR of meCpG ", breaks = 100)
dev.off()

## rule out IQR == 0 CpGs
idx_cpg_r <- is.na(cpg_val_iqr) | cpg_val_iqr == 0
cpg_r <- row.names(cpg_val)[idx_cpg_r]
idx_cor_cpg <- !is.na(match(cor_res$cpg, cpg_r))

## CTCF binding intensity 
ctcf_val <- read.table(gzfile("../Public_Data/ENCODE/ENCODE_ctcf_matrix.txt.gz"), header = T, row.names = 1)
ctcf_val_iqr <- apply(ctcf_val, 1,  IQR, na.rm = T)

pdf(file = "CTCF_intensity_IQR.pdf", width = 5, height = 3)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(ctcf_val_iqr, main = "", xlab = "IQR of CTCF intensity", breaks = 150)
dev.off()

# rule out IQR <1 0 CTCFs
idx_ctcf_r <- is.na(ctcf_val_iqr) | ctcf_val_iqr < 1
ctcf_r <- row.names(ctcf_val)[idx_ctcf_r]
idx_cor_ctcf <- !is.na(match(cor_res$ctcf_id, ctcf_r))

## filtering out meCpG-CTCP pairs based on samples cnt, CTCF and CpG IQR 
idx_r1 <- is.na(cor_res$pcc_2) | cor_res$numTested_o20 < 10 | idx_cor_cpg | idx_cor_ctcf
cor_res_o20 <- cor_res[!idx_r1, ]

ctcf_num <- length(unique(cor_res_o20$ctcf_id))
cpg_num <- length(unique(cor_res_o20$cpg))
dim(cor_res_o20)

pdf(file = "Number_of_encode_samples_for_CTCF_with_meCpG.pdf", width = 4, height = 3)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(cor_res_o20$numTested_o20, main ="with_meCpGs_inCTCF", xlab = "#ENOCDE tissue/cell-line")
dev.off()

pdf(file = "Number_of_CpGs_per_CTCF.pdf", width = 4, height = 3)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(table(cor_res_o20$ctcf_id), main = paste0("Number of CTCF: ", ctcf_num), xlab = "#CpGs", breaks = 40)
dev.off()
}

############################
## Multiple test corrections
cor_res_o20_out <- data.frame(cor_res_o20[, c(1:2,6:8)])        ## ussing pcc_2 and pval_2
qval_global <- p.adjust(cor_res_o20_out$pval_2, method = "BH")

##################################
## adding cpg and ctcf information 
{
  ## split CpG information
  ## 2 bases bed files cor_res_o20_sig
  cpg_info <- unlist(strsplit(cor_res_o20_out$cpg, ":"))
  chr = cpg_info[seq(1, length(cpg_info), 2)]
  
  ## 2 bases for bed files 
  cpg_start = as.numeric(cpg_info[seq(2, length(cpg_info), 2)]) - 1
  cpg_end = as.numeric(cpg_info[seq(2, length(cpg_info), 2)]) + 1 
  cpg_pos = as.numeric(cpg_info[seq(2, length(cpg_info), 2)])
  cpg_bed <- data.frame(chr, cpg_start, cpg_end, cpg_pos)
  
  ## adding the CTCF coordinates infomation 
  ctcf_bed <- read.table("../Public_Data/ENCODE/ENCODE_ctcf.bed")
  idx_ctcf <- match(cor_res_o20_out$ctcf_id, ctcf_bed$V4)
  
  ctcf_start = ctcf_bed$V2[idx_ctcf] 
  ctcf_end = ctcf_bed$V3[idx_ctcf]
  ctcf_mid = (ctcf_start + ctcf_end) / 2
  cpg2ctcf_mid = cpg_pos - ctcf_mid
  
  cpg_ctcf <- cbind(ctcf_start, ctcf_end, ctcf_mid, cpg2ctcf_mid)
}

cor_res_o20_out <- data.frame(cor_res_o20_out, qval_global, cpg_bed, cpg_ctcf)
#write.table(cor_res_o20_out, file = "cor_res_o20_all.txt", quote = F, row.names = F, col.names = T, sep = "\t")
}


####################################
## significant correlated meCpG-CTCF 
####################################
{
idx_sig <- abs(cor_res_o20_out$pcc_2) >= 0.5  & cor_res_o20_out$qval_global < 0.05
idx_sig_p <- idx_sig &  cor_res_o20_out$pcc_2 > 0.5
idx_sig_n <- idx_sig &  cor_res_o20_out$pcc_2 < 0.5

## plot pcc vs padj for all 
dat <- cor_res_o20_out
g <- ggplot(dat, aes(x = pcc_2, y = -log10(qval_global))) + geom_hex(bins = 50, show.legend = T)
g <- g + geom_hline(yintercept = 1, linetype = "dashed", color = "brown3")
g <- g + geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "brown3")
g <- g + labs(x = "cor(meCpG, CTCF)", y = "-log10(FDR)") + theme_classic()
ggsave("cor_res_o20_PCC_vs_padj_global_volcano.pdf", width = 4, height = 3)

## plot distance between meCpG to the center of CTCF
cor_dir <- rep("notSig", length(idx_sig))
cor_dir[idx_sig_p] <- "Pos"
cor_dir[idx_sig_n] <- "Neg"

dat <- data.frame(cor_res_o20_out, cor_dir)
g <- ggplot(dat, aes(x=cpg2ctcf_mid, color=cor_dir)) + geom_density() 
g <- g + xlab("Distance between meCpG and CTCF") + guides(color=guide_legend(title="cor(meCpG, CTCF)")) 
g <- g + theme_classic() + theme(legend.position="top")   
g <- g + scale_color_manual(values =  c("coral2", "gray", "cadetblue3"))
ggsave("DensityPlot_4_distance_between_CpG_and_CTCF_middle_point.pdf", width = 4, height = 4)

## K-S test 
ks.test(cor_res_o20_out$cpg2ctcf_mid[!idx_sig], cor_res_o20_out$cpg2ctcf_mid[idx_sig_n])
ks.test(cor_res_o20_out$cpg2ctcf_mid[!idx_sig], cor_res_o20_out$cpg2ctcf_mid[idx_sig_p])
ks.test(cor_res_o20_out$cpg2ctcf_mid[idx_sig_n], cor_res_o20_out$cpg2ctcf_mid[idx_sig_p])

## save the significant correlated pairs
cor_res_o20_sig <- cor_res_o20_out[idx_sig, ]
#write.table(cor_res_o20_sig, file = "cor_res_o20_sig_all.txt", quote = F, row.names = F, col.names = T, sep = "\t")

###############################################################
## filter out inconsistent correlated meCpG-CTCF pairs per CTCF
## number of meCpG per CTCF
pdf(file = "Number_of_sig_cor_CpGs_per_CTCF_with_meCpG_o20_all.pdf", width = 4, height = 3.5)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(table(cor_res_o20_sig$ctcf_id), main = "", xlab = "Number of meCpG per CTCF", breaks = 50)
dev.off()

## check whether PCC direction are consistent per CTCF
cor_res_o20_sig <- group_by(cor_res_o20_sig, ctcf_id) %>%  mutate(cor_cpg_ctcf_consis = mean(sign(pcc_2)))
idx_consist <- cor_res_o20_sig$cor_cpg_ctcf_consis == 1 | cor_res_o20_sig$cor_cpg_ctcf_consis == -1 
cor_res_o20_sig$cor_cpg_ctcf_consis[!idx_consist] = 0

length(cor_res_o20_sig$ctcf_id[!idx_consist])
length(unique(cor_res_o20_sig$ctcf_id[!idx_consist]))

## filter out non consistent ones
cor_res_o20_sig <- cor_res_o20_sig[idx_consist, ]
dim(cor_res_o20_sig)
length(unique(cor_res_o20_sig$ctcf_id))    ## number of CTCF

pdf(file = "Number_of_sig_cor_CpGs_per_CTCF_with_meCpG_o20_filtered.pdf", width = 4, height = 3)
par(mar = c(3, 3, 2, 1), mgp = c(2, 0.75, 0))
hist(table(cor_res_o20_sig$ctcf_id), main ="#sig_meCpGs_perCTCF", xlab = "#CpGs", breaks = 50)
dev.off()

# write.table(cor_res_o20_sig, file = "cor_res_o20_sig_filtered.txt", quote = F, row.names = F, col.names = T, sep = "\t")

### bed file for sig cor(CG, CTCF)
sig_cpg_bed <-  ungroup(cor_res_o20_sig) %>%  
                select(chr, cpg_start, cpg_end, cpg, ctcf_id) %>%  
                mutate(strand = ".", cor_cpg_ctcf = cor_res_o20_sig$pcc_2,
                                     cor_pval = cor_res_o20_sig$pval_2, 
                                     cor_fdr = cor_res_o20_sig$qval_global )

# write.table(sig_cpg_bed, file = "cor_res_o20_sig_filtered.bed", quote = F, row.names = F, col.names = T, sep = "\t")
}


############################################
### keep the most sig CpG per each CTCF peak 
#############################################
{
## grouping by CTCF
sig_cpg_bed_uniq <- group_by(sig_cpg_bed, ctcf_id) %>% 
                    filter(cor_pval == min(cor_pval)) %>% 
                    distinct(ctcf_id, .keep_all = TRUE)        ## remain only only CpG site even few with same p val

write.table(sig_cpg_bed_uniq, file = "cor_res_o20_sig_filtered_mostSig_cpg_perCTCF.bed", quote = F, row.names = F, col.names = T, sep = "\t")

idx_sig_n <-  sig_cpg_bed_uniq$cor_cpg_ctcf < -0.5
idx_sig_p <-  sig_cpg_bed_uniq$cor_cpg_ctcf > 0.5

dat <- sig_cpg_bed_uniq
g <- ggplot(dat, aes(x = cor_cpg_ctcf, y = -log10(cor_fdr))) + geom_hex(bins = 50, show.legend = T)
g <- g + geom_hline(yintercept = 1, linetype = "dashed", color = "brown3")
g <- g + geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "brown3")
g <- g + labs(x = "cor(meCpG, CTCF)", y = "-log10(FDR)") + theme_classic()
ggsave("cor_res_o20_PCC_vs_padj_global_mostSig_cpg_perCTCF_volcano.pdf", width = 4, height = 3)

#### adding corresponding information 
idx_mm <- match(sig_cpg_bed_uniq$cpg, cor_res_o20_sig$cpg)
cor_res_o20_sig_uniq <- cor_res_o20_sig[idx_mm, ]
write.table(cor_res_o20_sig_uniq, file = "cor_res_o20_sig_filtered_mostSig_cpg_perCTCF.txt", quote = F, row.names = F, col.names = T, sep = "\t")


##################################################################
## plot distance between cpg and middle point of CTCF
cor_dir <- rep("Pos", length(idx_sig_p))
cor_dir[idx_sig_n] <- "Neg"

dat <- data.frame(cor_res_o20_sig_uniq, cor_dir)
g <- ggplot(dat, aes(x=cpg2ctcf_mid, color=cor_dir)) + geom_density() 
g <- g + xlab("Distance between meCpG and CTCF") + guides(color=guide_legend(title="cor(meCpG, CTCF)")) 
g <- g + theme_classic() + theme(legend.position="top")  
ggsave("DensityPlot_4_distance_between_CpG_and_CTCF_middle_point_mostSig_cpg_perCTCF.pdf", width = 4, height = 4)

## K-S test 
ks.test(cor_res_o20_sig_uniq$cpg2ctcf_mid[idx_sig_n], cor_res_o20_sig_uniq$cpg2ctcf_mid[idx_sig_p])
}


##############################################
## comparison of mean meCpG and CTCF intensity
##############################################
{
idx_ctcf <- match(cor_res_o20_sig_uniq$ctcf_id, rownames(ctcf_val))
idx_cpg <- match(cor_res_o20_sig_uniq$cpg, rownames(cpg_val))

ctcf_mean <- rowMeans(ctcf_val[idx_ctcf, ], na.rm = T)
cpg_mean <- rowMeans(cpg_val[idx_cpg, ], na.rm = T)

group <- cor_dir 
group_col <- c("coral2", "cadetblue3")
cmp <- list(c("Neg", "Pos"))

dat <- data.frame(ctcf_mean, cpg_mean, group)

## for CTCF_mean
g1 <- ggviolin(data = dat, x = "group", y = "ctcf_mean", fill = "group", xlab = "Group", ylab = "Mean CTCF-binding intensity",
               palette = group_col, add = "boxplot", add.params = list(fill = "white"))
g1 <- g1 +  theme(legend.position = "none") # + stat_compare_means(comparisons = cmp)
ggsave("Mean_CTCF_signal_Cor_Neg_vs_Pos.pdf", width = 4, height = 4)


## for CpG_mean
g2 <- ggviolin(data = dat, x = "group", y = "cpg_mean", fill = "group", xlab = "Group", ylab = "Mean CpG methylation level",
               palette = group_col, add = "boxplot", add.params = list(fill = "white"))
g2 <- g2 + theme(legend.position = "none") # + stat_compare_means(comparisons = cmp)
ggsave("Mean_mCpG_signal_Cor_Neg_vs_Pos.pdf", width = 4, height = 4)
}
