##########################################################################
## Examine the relationship between memo-eQTLs and chormatin 3D structure 
## generating figures for Fig4 and Fig4S
rm(list = ls())
setwd("./Results")

library(GenomicRanges)
library(dplyr)
library(ggplot2)

#############
## memo-eQTLs
#############
res <- readRDS("./memo-eQTLs_meta_info_1711.RDS")

##########
## sig both
{
idx_sig_both <- res$sigGroup == "sigBoth"
res_both <- res[idx_sig_both & res$clus1.coef.svar.L * res$clus2.coef.svar.L < 0,  ]    ## removed one case
idx_g1 <- res_both$clus1.coef.svar.L < 0 &  res_both$clus2.coef.svar.L > 0
idx_g2 <- res_both$clus1.coef.svar.L > 0 &  res_both$clus2.coef.svar.L < 0

sigBoth_group <- vector()
sigBoth_group[idx_g1] <- "Group 1"
sigBoth_group[idx_g2] <- "Group 2"
res_both <- cbind(res_both, sigBoth_group)
#saveRDS(res_both, file = "memo-eQTLs_meta_info_sigBoth.rds")
dim(res_both)
}

## snp-gene range 
snp_gene_pos <- select(res, snp_pos, tss_pos) 
left <- apply(snp_gene_pos, 1, min)
right <- apply(snp_gene_pos, 1, max) 
snp_gene_range <- GRanges(res$chr, IRanges(left, right, names = rownames(res)))

cpg_loci <- GRanges(res$chr, IRanges(res$cpg_pos, res$cpg_pos + 1, names = rownames(res)))

## CTCF peaks range
ctcf <- read.table("../Public_Data/ENCODE/ENCODE_ctcf.bed")
idx_1 <- match(res$ctcf_id, ctcf$V4)
ctcf_loci <- GRanges(ctcf$V1[idx_1], IRanges(ctcf$V2[idx_1], ctcf$V3[idx_1], names = ctcf$V4[idx_1]))


##################################################
### check the overlapping with CTCF and/or H3K27ac 
### RWPE1 : normal prostate ChIA-PET
{
  IN1 <- read.table("../Public_Data/3D/LHR0001H.no_input_all_peaks.maxscore1000.narrowPeak")
  ctcf_1 <- GRanges(IN1$V1, IRanges(IN1$V2, IN1$V3))
 
  cpg_ctcf_idx <- findOverlaps(cpg_loci, ctcf_1)
  length(unique(from(cpg_ctcf_idx)))
  
  ctcf_ctcf_idx <- findOverlaps(ctcf_loci, ctcf_1)
  length(unique(from(ctcf_ctcf_idx)))
  
  
  ### read in CTCF ChIA-PET loop results 
  IN <- read.table("../Public_Data/3D/LHR0001H.e500.clusters.cis.BE3")
  ## requires at least 10 reads per CTCF loops
  #tmp <- IN[IN$V7 >= 10, ]    ## test different cut-offs, higher cut-off less meGpG in ctcf_loop
  
  tmp <- IN[IN$V7 >= 15, ]    ## similar to VcaP hi-c data 
  dim(tmp)
  
  ctcf_loop <- GRanges(tmp$V1, IRanges(tmp$V2, tmp$V6))
  ctcf_left <- GRanges(tmp$V1, IRanges(tmp$V2, tmp$V3))
  ctcf_right <- GRanges(tmp$V4, IRanges(tmp$V5, tmp$V6))
  
  ctcf_out <- as.data.frame(ctcf_loop)
  write.table(ctcf_out, file = "/Users/yong/OneDrive - UHN/Projects/me-eQTL/changhai/3_3D_loop/ChIA-PET/CTCF_ChIA_PET_loop_reads_gt_5.bed",
              row.names = F, col.names = F, quote = F, sep = "\t")
  
  ## ctcf_loop  length distribution 
  pdf("CTCF_loop_ranges_distribution_RWPE1.pdf", width = 5, height = 3)
  plot(density(log10(tmp$V6 - tmp$V2)), xlab = "log10loop_len", main = "CTCT_Loop_RWPE1")
  dev.off()
  
  ## snp_gene range overlappe with ctcf loop
  snp_gene_ctcf_loop <- findOverlaps(snp_gene_range, ctcf_loop)
  length(unique(from(snp_gene_ctcf_loop)))
  
  ctcf_loop_within  <- findOverlaps(ctcf_loop, snp_gene_range, type = "within" )
  
  tt <- paste0(length(unique(from(snp_gene_ctcf_loop ))), " of ", 
               length(snp_gene_range), " snp-cpg-gene ranges overlap with ",
               length(to(snp_gene_ctcf_loop)), " CTCF loops")
  
  st <- paste0(length(unique(from(ctcf_loop_within))), " of ", 
               length(unique(to(snp_gene_ctcf_loop))), 
               " Unique CTCF loops are within snp-cpg-gene ranges")
  
  ## snp-gene range only overlap with single CTCF loop
  sum(table(from(snp_gene_ctcf_loop )) == 1) / nrow(res)
  
  pdf("snp-cpg-gene ranges overlap with RWPE1 CTCF loops.pdf", width = 5, height = 3)
  hist(table(from(snp_gene_ctcf_loop )), breaks = 40, main = tt, cex.main = 0.75, sub = st,
       xlab = "#overlapped CTCF loops per snp-cpg-gene range" )
  dev.off()
  
  pdf("snp-cpg-gene ranges overlap with RWPE1 CTCF loops_label_adj.pdf", width = 4.5, height = 4)
  hist(table(from(snp_gene_ctcf_loop )), breaks = 40, main = "", col = "skyblue",
       xlab = "Number CTCF loops per eSNP-eCpG-eGene" )
  dev.off()
  
  ## in contrast
  pdf("snp-cpg-gene ranges overlap with RWPE1 CTCF loops per CTCF.pdf", width = 5, height = 3)
  hist(table(to(snp_gene_ctcf_loop )), breaks = 40, col = "skyblue",
       xlab = "#overlapped snp-cpg-gene range per CTCF loop", main = "" )
  dev.off()
  
  ##########################################################
  ## how many sigCpG sites are on the CTCF loop anchor sites
  cpg_ctcf_left <- findOverlaps(cpg_loci, ctcf_left)
  cpg_ctcf_right <- findOverlaps(cpg_loci, ctcf_right)
  length(unique(from(cpg_ctcf_left)))
  length(unique(from(cpg_ctcf_right)))
  
  ## CTCF loops within sigGpG
  ctcf_loop_cpg <- ctcf_loop[c(to(cpg_ctcf_left), to(cpg_ctcf_right)), ]
  names(ctcf_loop_cpg) <- c(names(cpg_loci)[from(cpg_ctcf_left)], names(cpg_loci)[from(cpg_ctcf_right)])
  
  length(ctcf_loop_cpg)                        ### all CTCF loop with memo-CpG sites
  length(unique(names(ctcf_loop_cpg)))         ## number of unique memo-eQTLs overlapped with CTCF loops
  
  ## same cpg site might overlap with multiple ctcf anchors !!!
  within_tmp <- findOverlaps(ctcf_loop_cpg, snp_gene_range, type = "within")
  
  ## making sure the ctcf_loop_cpg within the corresponding snp_gene_range
  idx_m <- names(ctcf_loop_cpg)[from(within_tmp)] == names(snp_gene_range)[to(within_tmp)]
  within_tmp <-  within_tmp[idx_m, ]
  
  length(unique(from(within_tmp)))
  within_name <- names(ctcf_loop_cpg)[from(within_tmp)]
  
  ####################
  ###### subset of res 
  ## selected unique ctcf_loop_cpg   
  tmp_cnt <- table(names(ctcf_loop_cpg))
  sum(tmp_cnt)
  length(unique(names(tmp_cnt )))
  
  pdf("RWPE1 CTCF loops shared anchor meGpG.pdf", width = 5, height = 3)
  hist(tmp_cnt, breaks = 20,
       xlab = "#CTCF loops share same meGPG", main = "" )
  dev.off()
  
  #idx_cnt <- tmp_cnt == 1  ## but the snp-gene range might still overlapped with other none ctcf_loop_cpg

  idx_cnt <- tmp_cnt == 1
  ctcf_loop_cpg_un <- names(tmp_cnt)[idx_cnt]
  
  idx <- match(rownames(res), ctcf_loop_cpg_un)   ## rm duplication
  res_s <- res[!is.na(idx), ]
  
  overlap <- rep("Cross", nrow(res_s))
  idx_with <- match(rownames(res_s), within_name)
  overlap[!is.na(idx_with)] <- "Within"
  table(overlap)
  
  dat <- data.frame(res_s, overlap)
  write.csv(dat, file = "memo-eQTLs_meta_info_with_CTCF_loop_RWPE1.csv")
  saveRDS(dat, file = "memo-eQTLs_meta_info_with_CTCF_loop_RWPE1.rds")
  
  
  g <- ggplot(dat, aes(x = -log10(clus1.modelpval.value), y = -log10(clus2.modelpval.value),
                       color = sigGroup, shape = overlap))
  g <- g + geom_point() + labs(x = "-log10(pval of relL)", y = "-log10(pval of relH)")
  g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray40")
  g <- g + geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "gray40")
  g <- g + scale_color_manual(values = c("coral4", "coral", "cadetblue", "gray"))
  g <- g + theme_classic()
  g <- ggsave("memo-eQTLs_meLow_vs_meHigh_with_CTCT_loops_RWPE1.pdf", heigh = 3.5, width = 4.5)
  
  ## summarise count 
  cnt <- dat %>%  
    group_by(sigGroup, overlap) %>% 
    summarise(n = n())
  
  ## compare high vs low 
  #M <- matrix(c(14, 64,11, 72), 2, 2)
  M <- matrix(c(cnt$n[3:6]), 2, 2)
  rownames(M) <- c("Cross", "Within")
  colnames(M) <- c("sigHigh", "sigLow")
  chisq.test(M)
}


#################################
## CTCF-loop in VCaP using HiChIP 
#################################

{
### load CTCF peaks in VCaP 
IN1 <- read.table("../Public_Data/3D/VCaP_CTCF_EtOH_hg38.bed")
ctcf_1 <- GRanges(IN1$V1, IRanges(IN1$V2, IN1$V3))

## Vehicle only                                    
## checking memo-CpG in ctcf peaks
cpg_ctcf_idx <- findOverlaps(cpg_loci, ctcf_1)
length(unique(from(cpg_ctcf_idx)))

### read in CTCF HiChIP loop results in VeH treatment
tmp <- read.table("../Public_Data/3D/vcap_ctcf_VeH_loop_hg38.bed")
ctcf_loop <- GRanges(tmp$V1, IRanges(tmp$V2, tmp$V3))

tmp <- read.table("../Public_Data/3D/vcap_ctcf_VeH_loop_left_anchor_hg38.bed")
ctcf_left <- GRanges(tmp$V1, IRanges(tmp$V2, tmp$V3))

tmp <- read.table("../Public_Data/3D/vcap_ctcf_VeH_loop_right_anchor_hg38.bed")
ctcf_right <- GRanges(tmp$V1, IRanges(tmp$V2, tmp$V3))

## snp_gene range overlappe with ctcf loop
snp_gene_ctcf_loop <- findOverlaps(snp_gene_range, ctcf_loop)
length(unique(from(snp_gene_ctcf_loop)))

ctcf_loop_within  <- findOverlaps(ctcf_loop, snp_gene_range, type = "within" )

tt <- paste0(length(unique(from(snp_gene_ctcf_loop ))), " of ", 
             length(snp_gene_range), " snp-cpg-gene ranges overlap with ",
             length(to(snp_gene_ctcf_loop)), " CTCF loops")

st <- paste0(length(from(ctcf_loop_within)), " of ", 
             length(to(snp_gene_ctcf_loop)), 
             " CTCF loops are within snp-cpg-gene ranges")

pdf("snp-cpg-gene ranges overlap with VCaP CTCF VeH loops.pdf", width = 5, height = 3)
hist(table(from(snp_gene_ctcf_loop )), breaks = 40, main = tt, cex.main = 0.75, sub = st,
     xlab = "#overlapped CTCF loops per snp-cpg-gene range" )
dev.off()


pdf("snp-cpg-gene ranges overlap with VCaP CTCF VeH loops_label_adj.pdf", width = 4, height = 4)
hist(table(from(snp_gene_ctcf_loop )), breaks = 40, main = "", col = "skyblue",
     xlab = "Number CTCF loops per eSNP-eCpG-eGene" )
dev.off()

## in contrast
pdf("snp-cpg-gene ranges overlap with VCaP CTCF VeH loops per CTCF.pdf", width = 5, height = 3)
hist(table(to(snp_gene_ctcf_loop )), breaks = 40,
     xlab = "#overlapped snp-cpg-gene range per CTCF loop", main = "" )
dev.off()


## how many sigCpG sites are on the CTCF loop anchor sites
cpg_ctcf_left <- findOverlaps(cpg_loci, ctcf_left)
cpg_ctcf_right <- findOverlaps(cpg_loci, ctcf_right)
length(from(cpg_ctcf_left))
length(from(cpg_ctcf_right))


## CTCF loops within snp-cpg-gene ranges 
ctcf_loop_cpg <- ctcf_loop[c(to(cpg_ctcf_left), to(cpg_ctcf_right)), ]
names(ctcf_loop_cpg) <- c(names(cpg_loci)[from(cpg_ctcf_left)], names(cpg_loci)[from(cpg_ctcf_right)])
## same cpg site might overlap with multiple ctcf anchors !!!
within_tmp <- findOverlaps(ctcf_loop_cpg, snp_gene_range, type = "within")
length(unique(from(within_tmp)))
within_name <- names(ctcf_loop_cpg)[unique(from(within_tmp))]


###### subset of res 
idx <- match(rownames(res), names(ctcf_loop_cpg))   ## rm duplications
res_s <- res[!is.na(idx), ]

overlap <- rep("Cross", nrow(res_s))
idx_with <- match(rownames(res_s), within_name)
overlap[!is.na(idx_with)] <- "Within"
table(overlap)

dat <- data.frame(res_s, overlap)
#saveRDS(dat, file = "memo-eQTLs_meta_info_with_CTCF_loop_VCaP.rds")
write.csv(dat, file = "memo-eQTLs_meta_info_with_CTCF_loop_VCaP.csv")


g <- ggplot(dat, aes(x = -log10(clus1.modelpval.value), y = -log10(clus2.modelpval.value),
                     color = sigGroup, shape = overlap))
g <- g + geom_point() + labs(x = "-log10(pval of relL)", y = "-log10(pval of relH)")
g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray40")
g <- g + geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "gray40")
g <- g + scale_color_manual(values = c("coral4", "coral", "cadetblue", "gray"))
g <- g + theme_classic()
g <- ggsave("memo-eQTLs_meLow_vs_meHigh_with_CTCT_loops.pdf", heigh = 3.5, width = 4.5)


## summarise count 
cnt <- dat %>%  
       group_by(sigGroup, overlap) %>% 
       summarise(n = n())

## compare high vs low 
M <- matrix(cnt$n[3:6], 2, 2)
rownames(M) <- c("Cross", "Within")
colnames(M) <- c("sigHigh", "sigLow")
chisq.test(M)
}
