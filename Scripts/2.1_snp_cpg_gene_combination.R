################################################################
## meCpG, SNP and gene expression filtering in 128 CPGEA samples
## identify all cpg-snp-gene combinations for memo-eQTL mapping
rm(list = ls())
setwd("./Results")

##############################################################
##       CPGEA samples' methylation data and filtering 
##############################################################
{
IN <- read.table("../Public_Data/CPGEA/CPGEA_samples.txt", as.is = T)
sample_id <- IN$V1
M <- length(sample_id)
  
raw <- read.table("../Public_Data/CPGEA/CPGEA_methylation.txt", as.is = T)
cpg_id <- paste(raw$V1, raw$V2, sep = ":")
N <- length(cpg_id)

me_beta <- matrix(0, N, M)                     ## cpg methylation ratio
me_tr <- matrix(0, N, M)                        ## cpg total reads
rownames(me_beta) <- rownames(me_tr) <- cpg_id
colnames(me_beta) <- colnames(me_tr) <- sample_id

## spiting samples, methylation levels
sample_sp <- strsplit(raw$V4, ",")
beta_sp <- strsplit(raw$V5, ",")
tr_sp <- strsplit(raw$V6, ",")

for(i in 1:N)
{
    idx <- match(sample_sp[[i]], sample_id)
    idx_s <- idx[!is.na(idx)]
    me_beta[i, idx_s] <- round(as.numeric(beta_sp[[i]])[!is.na(idx)], 2)         ## round 2 digits
    me_tr[i, idx_s] <- as.numeric(tr_sp[[i]])[!is.na(idx)]
}

## methylation median, mean and IQR
r_median <- apply(me_beta, 1, median)
r_mean <- rowMeans(me_beta)
r_iqr <- apply(me_beta, 1,  IQR)

pdf("Median_CpG_methy_in_128_Normal.pdf", width = 4, height = 3.5)
#plot(density(r_median), main = "Median_CpG_methy_in_128_Normal")
plot(density(r_median), main = "", xlab = "Median methylation level (Beta)", frame = FALSE)
dev.off()

### filter out CpG sites with less variable methylation levels
idx_me <- r_median > 0.25 & r_median < 0.75 
idx_iqr <- r_iqr > 0.1

### require significant correlated meCpG_CTCF in ENCODE      
cpg_ctcf <- read.table("../Data/cor_res_o20_sig_filtered_mostSig_cpg_perCTCF.txt", header = T)
idx_f <- match(rownames(me_beta), cpg_ctcf$cpg)

idx_k <- idx_me & idx_iqr & !is.na(idx_f)

pdf("Median_CpG_methy_in_128_Normal_filtered.pdf", width = 4, height = 4)
plot(density(r_median[idx_k]), main = "Median_CpG_methy_in_128_Normal_filtered")
dev.off()

write.table(me_beta[idx_k, ], file = "CPGEA_mdata_normal_filtered.txt", row.names = T, col.names = T, quote = F, sep = "\t")
#saveRDS(me_beta[idx_k, ], file = "CPGEA_normal_methy_beta_filtered.rds")
#saveRDS(me_tr[idx_k, ], file = "CPGEA_normal_methy_total_reads_filtered.rds")

## add cor (meCpG, CTCF) information 
me_f <- rownames(me_beta[idx_k, ])

idx_out <- match(me_f, cpg_ctcf$cpg)
cpg_ctcf_f <- cpg_ctcf[idx_out, ]

## cpg bed file 
cpg_bed <- cpg_ctcf_f %>%  select(chr, cpg_start, cpg_end, cpg, cpg_pos) %>%  mutate(strand = ".")
ctcf_bed <- cpg_ctcf_f %>% select(chr, ctcf_start, ctcf_end, ctcf_id, ctcf_mid) %>% mutate(strand = ".")

colnames(cpg_bed)[1] <- paste0("#", colnames(cpg_bed)[1])
colnames(ctcf_bed)[1] <- paste0("#", colnames(ctcf_bed)[1])

write.table(cpg_bed, file = "CPGEA_mdata_normal_filtered_CpG.bed", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(ctcf_bed, file = "CPGEA_mdata_normal_filtered_CTCF.bed", col.names = T, row.names = F, quote = F, sep = "\t")

}


##############################################################
##           CPGEA samples' expression data 
## filtering gene based on expression level and gene types
##############################################################
{
rpkm <- read.table("../Public_Data/CPGEA/CPGEA_FPKM.txt", header = T, as.is = T)       ## ENSG IDs might have same gene name
gene_gtf <- read.table("../Public_Data/GTF/GRCh38.84_Gene.gtf", as.is = T)
idx_s <- match(sample_id, colnames(rpkm))

## remove gene_name with at least two ensemble ID
idx_dup <- table(rpkm[, 2]) >= 2
gene_dup <- names(table(rpkm[, 2]))[idx_dup]
idx_rm <- match(rpkm[, 2], gene_dup)

rpkm_m <- rpkm[is.na(idx_rm), idx_s]
rownames(rpkm_m) <- rpkm[, 2][is.na(idx_rm)]

## checking gene types and keep protein_coding, lincRNA 
idx_m <- match(rownames(rpkm_m), gene_gtf$V4)
gene_type <- gene_gtf$V7[idx_m]
idx_tk <- gene_type == "protein_coding" | gene_type == "lincRNA" 
idx_tk[is.na(idx_tk)] <- FALSE

## protein_coding and lincRNA
rpkm_pe <- rpkm_m[idx_tk, ]
rpkm_pe_median  <- apply(rpkm_pe, 1, median)

pdf("Median_gene_expr_in_128_Normal.pdf", width = 4, height = 3.5)
plot(density(log2(rpkm_pe_median + 1)), xlab = "log2(Median FPKM + 1)", main = "", frame = FALSE)
dev.off()

## filtering based on the median expression level > 1 and gene type 
idx_re <- apply(rpkm_m, 1, median) > 1 & idx_tk 
rpkm_f <- rpkm_m[idx_re, ]

write.table(rpkm_f, file = "CPGEA_edata_normal_filtered.txt", quote = F, sep = "\t", row.names = T, col.names = T)

##### TSS site bed file for genes
idx_g <- match(rownames(rpkm_f), gene_gtf$V4) 
tss_bed <- gene_gtf[idx_g, c(1:4, 7, 6)]
idx_plus <- tss_bed$V6 == "+"

## plus strand
tss_bed[idx_plus, 3] <-  tss_bed[idx_plus, 2]
tss_bed[idx_plus, 2] <-  tss_bed[idx_plus, 2] - 1

## minus strand 
tss_bed[!idx_plus, 2] <-  tss_bed[!idx_plus, 3] - 1

colnames(tss_bed) <- c("#chr", "TSS_pos-1", "TSS_pos", "Gene", "Gene_type", "strand")
write.table(tss_bed, file = "CPGEA_edata_normal_filtered_TSS.bed", quote = F, sep = "\t", row.names = F, col.names = T)
}


##############################################################
##        CPGEA samples' genotype data and filtering 
##############################################################
{
## focusing on SNPs with rsID
gt_ori <- read.table("../Public_Data/CPGEA/CPGEA_genotype.txt", header = T, row.names = 1)
idx_s <- match(sample_id, colnames(gt_ori))

###  snp in PRAD ATAC-seq distall peaks 
if(FALSE){
atac <- read.table("../Public_Data/TCGA/ATAC-seq_PRAD_distal_peaks.txt", as.is = T)
atac_bed <- atac %>%  mutate(chr = paste0("chr", V1), ATAC_start = V2, ATAC_end = V3, ATAC_id = V4, ATAC_socre = V5, strand = ".") %>% select(chr:strand)
colnames(atac_bed)[1] <- "#chr"
write.table(atac_bed, file = "../Public_Data/TCGA/ATAC-seq_PRAD_distal_peaks.bed", row.names = F, col.names = T, quote = F, sep = "\t")
## intersectBed -a CPGEA_genotype.bed  -b ATAC-seq_PRAD_peakCalls.distal.bed  >  CPGEA_gdata_rsID_in_ATAC.bed
}

gt_bed_ori <- read.table("./CPGEA_SNP_in_ATAC.bed", as.is = T)

idx_na <- rowSums(is.na(gt_ori)) > 0  ## rm snp gt with NA 
idx_in <- match(rownames(gt_ori), gt_bed_ori$V4)     ## snp in ATAC
idx_gk <- !idx_na  & !is.na(idx_in) 

gt <- gt_ori[idx_gk, idx_s]

##########################
## MAF filterring: MAF > 5%  
D <- 2*ncol(gt)
alt_f <- rowSums(gt)/D 
ref_f <- 1 - alt_f
ff <- cbind(alt_f, ref_f)
maf <- apply(ff, 1, min)
gt <- gt[maf > 0.05, ]

###################
## filtered snp bed file
idx_bed <- match(rownames(gt), gt_bed_ori$V4)
gt_bed <- gt_bed_ori[idx_bed, ]

## rsID
tmp <- unlist(strsplit(gt_bed$V4, "_"))
rsid <- tmp[seq(1, length(tmp), 2)]
rownames(gt) <- rsid   ## gt rsID

write.table(gt, file = "CPGEA_gdata_rsID_filtered.txt", quote = F, sep = "\t", row.names = T, col.names = T)

## SNP pruning using the PLINK with R2 cutoff equals to 0.8
## and obtain pruned files :: CPGEA_gdata_rsID_filtered_LD_pruned_0.8.bed

## intersectBed 
## intersectBed -a ../Data/CPGEA_gdata_rsID_filtered_LD_0.8_pruned.bed  -b ../Data/ATAC-seq_PRAD_distal_peaks.bed -wa -wb >  CPGEA_gdata_rsID_filtered_LD_0.8_pruned_in_ATAC.bed

snp_atac <- read.table("./CPGEA_gdata_rsID_filtered_LD_0.8_pruned_in_ATAC.bed", as.is = T)
atac <- table(snp_atac$V12)
pdf("Number_of_filtered_SNPs_in_ATAC_peaks.pdf", width = 5, height = 4)
hist(atac, breaks = 16, xlab = "Number of SNPs in ATAC", 
     main = paste0("SNPs: ", nrow(snp_atac), "\t", "ATAC_peak:", length(atac)))
dev.off()

length(unique(snp_atac$V4))      ## number of SNPs
length(unique(snp_atac$V12))     ## number of ATAC-seqs
}


###########################################
## combination for snp-CpG-gene for testing
###########################################
{
## function for all combinations SNP-CpG-Gene combos in CPG centered 1Mbp region
  cpg_snp_gene <- function(file_tss_loci, file_snp_loci, file_cpg_loci, ExtDistance = 500000)
  {
    library(ChIPpeakAnno)
    #######################################################
    ## read in expression, genotype, methylation matrix ...
    {
      ## loci bed files without header
      tss_loci <- read.table(file_tss_loci, header = F, as.is = T)
      snp_loci <- read.table(file_snp_loci, header = F, as.is = T)
      cpg_loci <- read.table(file_cpg_loci, header = F, as.is = T)
    }
  
    #######################################################
    ##  CpG sites
    {
      ## Extending to CpG site centered ExtDistance*2 region
      cpg_loci_ext <- GRanges(seqnames = cpg_loci$V1,
                              ranges = IRanges(start = cpg_loci$V5 - ExtDistance,
                                               end = cpg_loci$V5 + ExtDistance,
                                               names = cpg_loci$V4))
      
      ## correct the pos < 0
      start(cpg_loci_ext)[start(cpg_loci_ext) < 0] <- 0
      
      cpg4test <- data.frame(cpg = cpg_loci$V4, chr = cpg_loci$V1, start = cpg_loci$V2, start = cpg_loci$V3)
 
    }
    
    #######################################################
    ##  SNPs and SNPs in CpG extented region
    {
      snp_loci_gr <- GRanges(snp_loci$V1, IRanges(snp_loci$V2, snp_loci$V3, names = snp_loci$V4))
      
      ## intersect SNPs with extented CpG region
      snp_in_cpg_ext <- subsetByOverlaps(snp_loci_gr, cpg_loci_ext)
      
      ## snp_in_cpg_ext for matrieqtl
      snp4test <- data.frame(snp = names(snp_in_cpg_ext), chr = seqnames(snp_in_cpg_ext), pos = end(snp_in_cpg_ext))
      
      ## Extending to SNP site centered "" region
      snp4test_ext <- ChIPpeakAnno::reCenterPeaks(snp_in_cpg_ext, ExtDistance*2)
      start(snp4test_ext)[start(snp4test_ext) < 0] <- 0
    }   
    
    
    #######################################################
    ##  SNP and genes combination
    {
      ## intersect Genes with extented SNP region
      tss_loci_gr <- GRanges(tss_loci$V1, IRanges(tss_loci$V2, tss_loci$V3, names = tss_loci$V4))
      
      ## snp list with all gene in extented region
      tss_in_snp_ext <- as.data.frame(findOverlaps(query = tss_loci_gr, subject = snp4test_ext))
      snps_genes_assoc <- split(tss_in_snp_ext$queryHits, as.character(tss_in_snp_ext$subjectHits))      ## genes in extended SNP region
      
      ## SNP rsIDs
      snp_idx <- as.numeric(names(snps_genes_assoc))
      names(snps_genes_assoc) <- names(snp4test_ext)[snp_idx]
      
      ## gene names
      snps_genes_assoc <- lapply(snps_genes_assoc, function(x) names(tss_loci_gr[x]))
    
      gene4test <- data.frame(gene_name = names(tss_loci_gr), chr = seqnames(tss_loci_gr), start = start( tss_loci_gr), end = end(tss_loci_gr))
      row.names(gene4test) <- gene4test$gene_name
    }
    
    
    #######################################################
    ## CpG and SNP combination
    {
      ## snps in CpG extended region
      snp_cpg_combo <- findOverlaps(query = snp_in_cpg_ext, subject = cpg_loci_ext)
      
      cpg_num <- length(unique(to(snp_cpg_combo)))     ## number of CpG sits with SNP in extended region
      cpg_idx <- to(snp_cpg_combo)
      
      cpgs_snps_assoc <- list()
      
      for (i in 1:cpg_num){
        j <- unique(cpg_idx)[i]
        cpgs_snps_assoc[[i]] <- names(snp_in_cpg_ext[from(snp_cpg_combo[cpg_idx == j])])
        names(cpgs_snps_assoc)[i] <- names(cpg_loci_ext)[j]
        }
      
    }
    
    ############################
    ## all possible combinations
    combo <- data.frame()
    L1 <- length(cpgs_snps_assoc)
    for (i in 1: L1)
    {
        ## matching snps
        idx_s <- match(cpgs_snps_assoc[[i]], names(snps_genes_assoc))
        
        if(sum(is.na(idx_s)) != length(idx_s))
        {
        snp_gene_tmp <- snps_genes_assoc[idx_s[!is.na(idx_s)]]
        snp_gene_v <- setNames(unlist(snp_gene_tmp, use.names=F), 
                               rep(names(snp_gene_tmp),
                               lengths(snp_gene_tmp)))
        cpg <-  rep(names(cpgs_snps_assoc)[i], length (snp_gene_v))
        snp <-  names(snp_gene_v)
        gene <- snp_gene_v
      
        combo_tmp <- data.frame(cpg, snp, gene)
        combo <- rbind(combo, combo_tmp)
        }
    }
    
    ## requiring the cpg in the middle of snp and gene
    ## add pos and distans info to combo
    idx_m <- match (combo$cpg, cpg4test$cpg)
    cpg_pos <- cpg4test$start[idx_m] + 1
    
    idx_m <- match (combo$snp, snp4test$snp)
    snp_pos <- snp4test$pos[idx_m]
    
    idx_m <- match (combo$gene, gene4test$gene_name)
    gene_pos <- gene4test$start[idx_m]
    
    snp2gene <- snp_pos - gene_pos
    cpg2snp <- cpg_pos - snp_pos
    cpg2gene <- cpg_pos - gene_pos
    idx_ss <-  cpg2snp * cpg2gene < 0
    
    combo <- cbind(combo, cpg_pos, snp_pos, gene_pos, snp2gene, cpg2snp, cpg2gene)
    ## range(snp2gene) 
    combo_ss <- combo[idx_ss, ]
    
    return(combo_ss)
  
  }


## run the function 
  {
    file_tss_loci <- "./CPGEA_edata_normal_filtered_TSS.bed"
    file_snp_loci <- "./CPGEA_gdata_rsID_filtered_LD_0.8_pruned_in_ATAC.bed"
    file_cpg_loci <- "./CPGEA_mdata_normal_filtered_CpG.bed"
    
    comb <- cpg_snp_gene(file_tss_loci = file_tss_loci,
                         file_snp_loci = file_snp_loci,
                         file_cpg_loci = file_cpg_loci)
    write.table(comb, file = "cpg_snp_gene_combination_for_memo-eQTL_mapping.txt",
                sep = "\t", row.names = F, col.names = F, quote = F)
  }
}

