rm(list = ls())
setwd("./Results")


library(dplyr)
library(stringr)
library(ggplot2)
## tested snp gene combinations
memo <- readRDS("./all_tested_snp_cpg_gene_cmb_4_memo-eqtl.RDS")

## M1 only and remove redundant SNP-Gene pairs due modulated by different CpG sites 
memo_m1 <- distinct(memo[, c(2:3, 7:12)])
m1_genes <- unique(memo_m1$evar)
m1_snps <- unique(memo_m1$svar)

####################################
## gene expression for GTEx prostate
tpm <- read.table(gzfile("../Public_Data/GTEx/gene_tpm_2017-06-05_v8_prostate.gct.gz"), 
                  skip = 2, header = T)
# gene name
idx_g <- match(m1_genes, tpm$Description)
tpm <- tpm[idx_g[!is.na(idx_g)], -c(1:3)]
rownames(tpm) <- m1_genes[!is.na(idx_g)]

# sample ids
n_split <- str_split_fixed(colnames(tpm), "\\.", 3)
n_short <- paste(n_split[, 1], n_split[, 2], sep = "-")
colnames(tpm) <- n_short

########################################################
## genotype for the SNPs in ATAC peaks for GTEx prostate
snp <- readRDS("../Public_Data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_Prostate_PASS.recode.GT.FORMAT_012_CPGEA_in_ATAC.RDS")

idx_snp <- match(m1_snps, rownames(snp))

## samples with both tpm and snp
idx_s <- match(colnames(snp), colnames(tpm))
snp_s <- as.data.frame(snp[idx_snp[!is.na(idx_snp)], !is.na(idx_s)])

tpm_s <- as.data.frame(tpm[, idx_s[!is.na(idx_s)]])

idx_gx <- match(memo_m1$evar, rownames(tpm_s)) 
idx_snpx <- match(memo_m1$svar, rownames(snp_s)) 

memo_m1_s <- memo_m1[!is.na(idx_gx) & !is.na(idx_snpx), ] 

dim(snp_s)
dim(tpm_s)
length(unique(memo_m1_s$evar))
length(unique(memo_m1_s$svar))

###########################################
## confirm eQTL fitting model
###########################################

if(FALSE){
## M1 validation 
memo_expr <- read.table("./CPGEA_edata_normal_filtered.txt", header = T)
memo_snp <- read.table("./CPGEA_gdata_rsID_filtered.txt", header = T) 
sum(colnames(memo_expr) == colnames(memo_snp))
memo_m1_s[1, ]
evar <- t(memo_expr[which(rownames(memo_expr) == memo_m1_s$evar[1]), ])
svar <- t(memo_snp[which(rownames(memo_snp) == memo_m1_s$svar[1]), ])

## without scaling 
df <- data.frame(evar, svar)
lm(evar ~ svar, df)

## with scaling 
evar_s <- apply(evar, 2, RNOmni::RankNorm)
svar_s <- svar / 2
df_s <- data.frame(evar_s, svar_s)
res <- summary(lm(evar_s ~ svar_s, df_s))
beta <- res$coefficients[2, 1]
pval <- res$coefficients[2, 4]
r2 <-   res$adj.r.squared
}


#####################################################
## perform eQTL by subsampling 128 samples
## for all SNP-Gene combinations tested for memo-eQTL
####################################################
if(FALSE){
k = 1
N <- nrow(memo_m1_s)
res_sampling <- list()

while( k <= 10)
{
set.seed(k)      ## ensure the sampling is reproducible 
  
idx_ss <- sample(1:ncol(snp_s), 128)
## scale snp to the range [0, 1]
snp_ss <- snp_s[, idx_ss] / 2
## apply RankNorm for gene expression 
tpm_ss <- as.data.frame(t(apply(tpm_s[, idx_ss], 1, RNOmni::RankNorm)))

beta <- pval <- rsq <- c()

for(i in 1:N){
  if(i %% 1000 == 0){
    print(paste0("Subsampling ", k, ": testing ", i, "of ", N, "SNP-Gene pairs ..."))
  }
  
  ## for M1 only 
  evar <- t(tpm_ss[which(rownames(tpm_ss) == memo_m1_s$evar[i]), ])
  svar <- t(snp_ss[which(rownames(snp_ss) == memo_m1_s$svar[i]), ])
  df <- data.frame(evar, svar)
  res <- summary(lm(evar ~ svar, df))
    if(nrow(res$coefficients) == 1){
      beta <- c(beta, NA)
      pval <- c(pval, NA)
      rsq  <- c(rsq,  NA)
    } else {
      beta <- c(beta, res$coefficients[2, 1])
      pval <- c(pval, res$coefficients[2, 4])
    rsq  <- c(rsq,  res$adj.r.squared)
    }
  }
res_sampling[[k]] <- data.frame(beta, pval, rsq)
k = k + 1
}

saveRDS(res_sampling, "GTEx_subsampling_128_eQTLs.RDS")
}

#####################################
## permutation results visualization
{
res_sampling <- readRDS("GTEx_subsampling_128_eQTLs.RDS")

###############################################
## check overlapped significant eQTL (p < 0.05)
idx_ref <- memo_m1_s$m1_svar_pval < 0.05

overlap_cnt <- c(sum(idx_ref), 0)

group <- c("Ref", "Ref", rep(paste0("Sampling_", 1:10), each = 2))
sub_group <- rep(c("Overlapped_with_Ref", "non-overlapped"), 11)

for(i in 1:10){
  idx_sig <- res_sampling[[i]][, 2] < 0.05
  idx_ov <- idx_sig & idx_ref 
  overlap_cnt <- c(overlap_cnt, sum(idx_ov, na.rm = T), sum(idx_sig, na.rm = T) - sum(idx_ov, na.rm = T))
}

dat <- data.frame(overlap_cnt, sub_group, group)

g <- ggplot(dat, aes(fill = sub_group, y = overlap_cnt, x = group, label = overlap_cnt))  
g <- g + geom_bar(position='stack', stat='identity')
g <- g + geom_text(size = 5, position = position_stack(vjust = 0.5))
g <- g + theme_classic()
ggsave("sig_eQTL_overlapped_with_GTEx_subsampling_128_eQTLs.pdf", width = 12, height = 4)

#########################################
## for non-significant pairs overlapping
idx_ref <- memo_m1_s$m1_svar_pval >= 0.05

overlap_cnt <- c(sum(idx_ref), 0)

group <- c("Ref", "Ref", rep(paste0("Sampling_", 1:10), each = 2))
sub_group <- rep(c("Overlapped_with_Ref", "non-overlapped"), 11)

for(i in 1:10){
  idx_non_sig <- res_sampling[[i]][, 2] >= 0.05
  idx_ov <- idx_non_sig  & idx_ref 
  overlap_cnt <- c(overlap_cnt, sum(idx_ov, na.rm = T), sum(idx_non_sig, na.rm = T) - sum(idx_ov, na.rm = T))
}

dat <- data.frame(overlap_cnt, sub_group, group)

g <- ggplot(dat, aes(fill = sub_group, y = overlap_cnt, x = group, label = overlap_cnt))  
g <- g + geom_bar(position='stack', stat='identity')
g <- g + geom_text(size = 5, position = position_stack(vjust = 0.5))
g <- g + theme_classic()
ggsave("non-sig_eQTL_overlapped_with_GTEx_subsampling_128_eQTLs.pdf", width = 12, height = 4)


## gray 
g <- ggplot(dat, aes(fill = sub_group, y = overlap_cnt, x = group, label = overlap_cnt))  
g <- g + geom_bar(position='stack', stat='identity')
g <- g + geom_text(size = 5, position = position_stack(vjust = 0.5))
g <- g + scale_fill_manual(values = c("gray40", "gray80"))
g <- g + theme_classic() + theme(legend.position = "top") 
ggsave("non-sig_eQTL_overlapped_with_GTEx_subsampling_128_eQTLs_pub.pdf", width = 9, height = 5)


#########################
## for sig-memo-eQTL only
{
sig_memo <- readRDS("./sig_memo-eQTLs.RDS")

## unique snp_gene pair
sig_pairs <- distinct(sig_memo[, 2:3])
idx_sig_m <- match(paste0(sig_pairs$evar, sig_pairs$svar), paste0(memo_m1_s$evar, memo_m1_s$svar), nomatch = 0)

memo_m1_s_m <- memo_m1_s[idx_sig_m, ]

## for sig
idx_ref <- memo_m1_s_m$m1_svar_pval < 0.05

overlap_cnt <- c(sum(idx_ref), 0)

group <- c("Ref", "Ref", rep(paste0("Sampling_", 1:10), each = 2))
sub_group <- rep(c("Overlapped_with_Ref", "non-overlapped"), 11)

for(i in 1:10){
  idx_non_sig <- res_sampling[[i]][idx_sig_m, 2] < 0.05
  idx_ov <- idx_non_sig  & idx_ref 
  overlap_cnt <- c(overlap_cnt, sum(idx_ov, na.rm = T), sum(idx_non_sig, na.rm = T) - sum(idx_ov, na.rm = T))
}

dat <- data.frame(overlap_cnt, sub_group, group)

g <- ggplot(dat, aes(fill = sub_group, y = overlap_cnt, x = group, label = overlap_cnt))  
g <- g + geom_bar(position='stack', stat='identity')
g <- g + geom_text(size = 5, position = position_stack(vjust = 0.5))
g <- g + theme_classic()
ggsave("sig_eQTL_for_memo-eQTLs_overlapped_with_GTEx_subsampling_128_eQTLs.pdf", width = 12, height = 4)

## for nog-sig
idx_ref <- memo_m1_s_m$m1_svar_pval >= 0.05

overlap_cnt <- c(sum(idx_ref), 0)

group <- c("Ref", "Ref", rep(paste0("Sampling_", 1:10), each = 2))
sub_group <- rep(c("Overlapped_with_Ref", "non-overlapped"), 11)

for(i in 1:10){
  idx_non_sig <- res_sampling[[i]][idx_sig_m, 2] >= 0.05
  idx_ov <- idx_non_sig  & idx_ref 
  overlap_cnt <- c(overlap_cnt, sum(idx_ov, na.rm = T), sum(idx_non_sig, na.rm = T) - sum(idx_ov, na.rm = T))
}

dat <- data.frame(overlap_cnt, sub_group, group)

g <- ggplot(dat, aes(fill = sub_group, y = overlap_cnt, x = group, label = overlap_cnt))  
g <- g + geom_bar(position='stack', stat='identity')
g <- g + geom_text(size = 5, position = position_stack(vjust = 0.5))
g <- g + theme_classic()
ggsave("non-sig_eQTL_for_memo-eQTLs_overlapped_with_GTEx_subsampling_128_eQTLs.pdf", width = 12, height = 4)

}



}

##########################################################################
## perform sig memo-eQTL by subsampling to same subgroups samples (relH/L)
## to test whether the CpG modulation effects are by chance
## random sampling 128 samples >> random partition 128 sample in two group based sig-memo splition 
#########################################################################
if(FALSE){
sig_memo <- readRDS("/sig_memo-eQTLs.RDS")

## remove sig_memo without available gene expression and SNP in GTEX
idx_1 <- match(sig_memo$evar, row.names(tpm_s))
idx_2 <- match(sig_memo$svar, row.names(snp_s))

## M1 and subgroups based the CpG relH/L 
sig_memo_s <- distinct(sig_memo[!is.na(idx_1) & !is.na(idx_2), c(2:4, 7:12, 38:51)])
dim(sig_memo_s)    ## 1016 examined!!


## CpG methyaltion level cut off 
hl_cut <- as.numeric(str_split_fixed(sig_memo_s$clus2label, "_", 4)[, 2])

## methylation levels
meth <- read.table("./CPGEA_mdata_normal_filtered.txt", header = T)
idx_m <- match(sig_memo_s$mvar,rownames(meth))
meth_m <- meth[idx_m, ]
n_h <- rowSums(meth_m >= hl_cut)

###############################################
## number of samples in CpG relative high group 

set.seed(11)      ## ensure the sampling is reproducible 
idx_ss <- sample(1:ncol(snp_s), 128)
## scale snp to the range [0, 1]
snp_ss <- snp_s[, idx_ss] / 2
## apply RankNorm for gene expression 
tpm_ss <- as.data.frame(t(apply(tpm_s[, idx_ss], 1, RNOmni::RankNorm)))


N <- nrow(sig_memo_s)
res_subgroup <- list()

  for(i in 1:N){
    
    print(paste0("Testing ", i, " of ", N, " sig memo-eQTLs ..."))
    
    ## for M1 only 
    beta <- pval <- rsq <- c()
    evar <- t(tpm_ss[which(rownames(tpm_ss) == sig_memo_s$evar[i]), ])
    svar <- t(snp_ss[which(rownames(snp_ss) == sig_memo_s$svar[i]), ])
    df <- data.frame(evar, svar)
    
    res <- summary(lm(evar ~ svar, df))
    if(nrow(res$coefficients) == 1){
      beta <- NA
      pval <- NA
      rsq  <- NA
    } else {
      beta <- res$coefficients[2, 1]
      pval <- res$coefficients[2, 4]
      rsq  <- res$adj.r.squared
    }
    m1 <- data.frame(sig_memo_s[i, ], beta , pval ,  rsq)
    
    ## for randomly partitions 
    k = 1
    beta_h <- pval_h <- rsq_h <-  beta_l <- pval_l <- rsq_l <- c()
    while (k <= 1000){
    set.seed(11 + k)
    idx_sh <- sample(1:ncol(snp_ss), n_h[i])
    idx_sh
    
    df_h <- df[idx_sh, ]
    df_l <- df[-idx_sh, ]
    
    ## pseudo high
    res_h <- summary(lm(df_h[, 1] ~ df_h[, 2], df_h))
    if(nrow(res_h$coefficients) == 1){
      beta_h <- c(beta_h, NA)
      pval_h <- c(pval_h, NA)
      rsq_h  <- c(rsq_h,  NA)
    } else {
      beta_h <- c(beta_h, res_h$coefficients[2, 1])
      pval_h <- c(pval_h, res_h$coefficients[2, 4])
      rsq_h  <- c(rsq_h,  res_h$adj.r.squared)
    }
    
    ## pseudo low
    res_l <- summary(lm(df_l[, 1] ~ df_l[, 2], df_l))
    if(nrow(res_l$coefficients) == 1){
      beta_l <- c(beta_l, NA)
      pval_l <- c(pval_l, NA)
      rsq_l  <- c(rsq_l,  NA)
    } else {
      beta_l <- c(beta_l, res_l$coefficients[2, 1])
      pval_l <- c(pval_l, res_l$coefficients[2, 4])
      rsq_l  <- c(rsq_l,  res_l$adj.r.squared)
    }
    k = k + 1
    }
    
    subgroup <- data.frame( beta_h , pval_h ,  rsq_h , beta_l , pval_l ,  rsq_l)
    
    res_subgroup[[i]] <- list(m1, subgroup)
    
  
    } 
  
saveRDS(res_subgroup, "GTEx_subsampling_128_and_HL_partitiion_memo-eQTLs.RDS")

}

#####################################
## permutation results visualization
{
res_subgroup <- readRDS("GTEx_subsampling_128_and_HL_partitiion_memo-eQTLs.RDS")
## for memo-eQTLs
H_pval <- L_pval <- c()
## random partition.
H_pval_1000 <- L_pval_1000 <- data.frame(rep(0, 1000))
for (i in 1:length(res_subgroup)){
  H_pval <- c(H_pval, res_subgroup[[i]][[1]][1, 21])     ## clus2.pval
  L_pval <- c(L_pval, res_subgroup[[i]][[1]][1, 14])     ## clus1.pval
  H_pval_1000[, i] <- res_subgroup[[i]][[2]][, 2]
  L_pval_1000[, i] <- res_subgroup[[i]][[2]][, 5]
}

#########
##sigHigh
{
H_pval_1001z <- scale(rbind(H_pval, H_pval_1000))
idx_so <- H_pval < 0.05 & L_pval >= 0.05
so_z <- H_pval_1001z[, idx_so]
## remove all NA columns
idx_k <- !(colSums(is.na(so_z)) == nrow(so_z))
so_z <- so_z[, idx_k]

k = 1
ref <- data.frame()
while(k <= 1000){
ref <- rbind(ref, so_z[1, ])
k = k + 1
}
sig_cnt <- colSums(ref < so_z[-1, ], na.rm = T)  
all_cnt <- colSums(!is.na(so_z[-1, ]))

idx_sig <- sig_cnt/all_cnt > 0.95 

## sort by color and ref_pz
so_zr <- so_z[, idx_sig]
idx_zr <- order(so_zr[1, ])
so_zr <- so_zr[, idx_zr]

so_zb <- so_z[, !idx_sig]
idx_zb <- order(so_zb[1, ])
so_zb <- so_zb[, idx_zb]

so_zs <- cbind(so_zr, so_zb)
colnames(so_zs) <- 1:ncol(so_zs)
col_sig <- rep("blue", ncol(so_zs))
col_sig[1:ncol(so_zr)] <- "red"

## ref_point modifcation
so_zsp <- so_zs[1, ]
idx_pm <- so_zsp < -3
so_zsp[idx_pm] <- -3
pch_sig <- rep(16, ncol(so_zs))
pch_sig[idx_pm] <- 18

pdf("Permutation test for sgHigh memo-eQTLs.pdf", width = 12, height = 4)
par(mar=c(4,4,1,1))
x_lab <- paste0("sigHigh memo-eQTLs: ", ncol(so_zs), "; Permutation tests: Significant(red)=", ncol(so_zr),
               "; non-significant(blue)=", ncol(so_zb))
boxplot(so_zs, outline = FALSE, xaxt = "n", xlab = x_lab, ylab = "Z score of p values")
points(1:ncol(so_zs), so_zsp, pch = pch_sig, cex = 0.35, col = col_sig)
dev.off()
}

#########
##sigLow
{
  L_pval_1001z <- scale(rbind(L_pval, L_pval_1000))
  idx_so <- L_pval < 0.05 & H_pval >= 0.05
  so_z <- L_pval_1001z[, idx_so]
  ## remove all NA columns
  idx_k <- !(colSums(is.na(so_z)) == nrow(so_z))
  so_z <- so_z[, idx_k]
  
  k = 1
  ref <- data.frame()
  while(k <= 1000){
    ref <- rbind(ref, so_z[1, ])
    k = k + 1
  }
  sig_cnt <- colSums(ref < so_z[-1, ], na.rm = T)  
  all_cnt <- colSums(!is.na(so_z[-1, ]))
  
  idx_sig <- sig_cnt/all_cnt > 0.95 
  
  ## sort by color and ref_pz
  so_zr <- so_z[, idx_sig]
  idx_zr <- order(so_zr[1, ])
  so_zr <- so_zr[, idx_zr]
  
  so_zb <- so_z[, !idx_sig]
  idx_zb <- order(so_zb[1, ])
  so_zb <- so_zb[, idx_zb]
  
  so_zs <- cbind(so_zr, so_zb)
  colnames(so_zs) <- 1:ncol(so_zs)
  col_sig <- rep("blue", ncol(so_zs))
  col_sig[1:ncol(so_zr)] <- "red"
  
  ## ref_point modifcation
  so_zsp <- so_zs[1, ]
  idx_pm <- so_zsp < -3.5
  so_zsp[idx_pm] <- -3.5
  pch_sig <- rep(16, ncol(so_zs))
  pch_sig[idx_pm] <- 18
  
  pdf("Permutation test for sgLow memo-eQTLs.pdf", width = 12, height = 4)
  par(mar=c(4,4,1,1))
  x_lab <- paste0("sigLow memo-eQTLs: ", ncol(so_zs), "; Permutation tests: Significant(red)=", ncol(so_zr),
                  "; non-significant(blue)=", ncol(so_zb))
  boxplot(so_zs, outline = FALSE, xaxt = "n", xlab = x_lab, ylab = "Z score of p values")
  points(1:ncol(so_zs), so_zsp, pch = pch_sig, cex = 0.35, col = col_sig)
  dev.off()
}

}
