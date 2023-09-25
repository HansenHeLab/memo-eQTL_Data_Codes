###################################################################
## Perform correlation test for all possible meCpG-CTCF pairs
## Shell script/command and R script as below were used on cluster
###################################################################


######################################################################
## all meCpG-CTCF pairs were split to 100000 lines for job submission 

# split -l 100000  ./Public_Data/ENCODE/meCpG_CTCF_pairs.txt
# for i in splits/*; do sbatch --export=infile=$i run_cortest_submitter.sh ; done

###########################
## run_cortest_submitter.sh 
{
#!/bin/bash
#SBATCH -p all
#SBATCH -e ./slurm-%j.err
#SBATCH -o ./slurm-%j.out
#SBATCH -t 10:00:00
#SBATCH --mem=8192M

# cd ./Public_Data/ENCODE/
  
# module load R/3.3.0
# Rscript run_cortest.R $infile 
}

##################################
## run_cortest.R 
## apply spearman correlation test 
{
  args = commandArgs(trailingOnly=TRUE)
  infile <- read.table(args[1], header = T, sep = "\t")
  me <- read.table("encode_methylation_matrix.txt", header = T, row.names = 1, sep = "\t")
  ct <- read.table("encode_ctcf_matrix.txt", header = T, row.names = 1, sep = "\t")
  tr <- read.table("encode_totalRead_matrix.txt", header = T, row.names = 1, sep = "\t")
  res <- data.frame()
  for ( i in 1:nrow(infile)){
    meid <- as.character(infile[i,1])
    cid <- as.character(infile[i,2])
    mevals <- as.numeric(me[meid,])
    cvals <- as.numeric(ct[cid,])
    e1 <- NA
    p1 <- NA
    e2 <- NA
    p2 <- NA
    if ( length(!is.na(mevals) & !is.na(cvals)) >= 10){
      ptest <- try(cor.test(mevals[!is.na(mevals) & !is.na(cvals)],cvals[!is.na(mevals) & !is.na(cvals)], method = "spearman"))
      if (class(ptest) == "htest"){
        e1 <- ptest$estimate
        p1 <- ptest$p.value
      }
      o20 <- which(!is.na(mevals) & !is.na(cvals) & tr[meid,] >=20)
      if (length(o20) >= 10){
        ptest <- try(cor.test(mevals[o20], cvals[o20]))
        if (class(ptest) == "htest"){
          e2 <- ptest$estimate
          p2 <- ptest$p.value
        }
      }
    }
    res <- rbind(res, data.frame(cid,meid,length(mevals[!is.na(mevals) & !is.na(cvals)]), e1,p1,length(o20),e2,p2))
  }
  write.table(res, file = paste0("./out_",args[1], ".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
}

#############################
## splited ruslt were mereged
# cat x*txt > ./Results/cor_meCpG_CTCF_pearson.txt 
# cat x*txt > ./Results/cor_meCpG_CTCF_spearman.txt 

