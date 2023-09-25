
rm(list = ls())
setwd("./Results")
## setwd("/Users/yong/OneDrive - UHN/Projects/me-eQTL/changhai/memo-eQTL_Data_Codes/Results")

library(readr)        ## read_table2
library(data.table)
library(plyr)
library(future.apply)
library(future)
library(dplyr)
library(mosaic)

# library(distributions3)

####################
#### Functions ####
###################
{
  # Data Preparation Function
  cohortdf = function(fileloc=".", pattern = "filtered.txt"){
    # Load the text files: eqtl, gene and meth
    filelist = list.files(fileloc, pattern="*filtered.txt", full.names = TRUE)
    tarfilelist = filelist[grep(pattern, filelist)]
    
    dfcreate = function(mainfile){
      df = data.frame(mainfile, stringsAsFactors = F)
      names(df) = df[1,]
      df = df[-1,]
      df[,-1] =apply(df[,-1],2,as.numeric)
      return(df)
    }
    ## gene
    filename = tarfilelist[grep("edata",tarfilelist)]
    mainfile <- read.table(filename, header = F, fill = T)
    gene = dfcreate(mainfile)
    # str(gene)
    
    ## snp
    filename = tarfilelist[grep("gdata",tarfilelist)]
    mainfile <- read.table(filename, header = F, fill = T)
    snp = dfcreate(mainfile)
    
    ## meth
    filename = tarfilelist[grep("mdata",tarfilelist)]
    mainfile <- read.table(filename, header = F, fill = T)
    meth = dfcreate(mainfile)
    # str(meth)
    
    df = list(gene, snp, meth)
    return(df)
    return(df)
  }
  
  ## Genomic dataset preparation
  genedataprep = function(fileurl){
    
    # Get the names of the columns provided by Yong
    testfile <- read_table2(fileurl)
    testnames = data.frame(testfile, stringsAsFactors = F)
    for(i in 1:3){
      a= testnames[,i]
      if(all(grepl("chr", a))){names(testnames)[i] = "mvar"}
      else if(all(grepl("rs", a))){names(testnames)[i] = "svar"}
      else{names(testnames)[i] = "evar"}
    }
    testnames = testnames[,c("evar", "svar", "mvar")]
    df = testnames[!duplicated(testnames),]
    
    return(df)
  }
  
  ## Prepare the samples for each gene, snp and meth
  dataprep = function(maindf = Prefile, testdf){
    sampledf = lapply(maindf, function(x){
      y = data.frame(t(x[,-1]))
      names(y) = x[,1]; #str(y);
      y})
    sampledf = lapply(1:3, function(x){ y = sampledf[[x]];
    y = y[,names(y) %in% unique(testdf[,x])]
    # str(y)
    y
    })
    return(sampledf)
  }
  
  
  ## Remove names not found in the main sample
  cleangenedata = function(maindf = Prefile, testdf, savename = "test"){
    
    sampledf = dataprep(maindf = maindf, testdf = testdf)
    
    # find and save names not found in main sample:
    evarnf = setdiff(testdf$evar , names(sampledf[[1]]))
    nfd = data.frame(gene = evarnf)
    filename = paste("genenf",savename,".csv", sep="")
    write.csv(nfd, filename)
    
    svarnf = setdiff(testdf$svar, names(sampledf[[2]]))
    nfd = data.frame(snp = svarnf)
    filename = paste("snpnf",savename,".csv", sep="")
    write.csv(nfd, filename)
    
    mvarnf = setdiff(testdf$mvar, names(sampledf[[3]]))
    nfd = data.frame(meth = mvarnf)
    filename = paste("methnf",savename,".csv", sep="")
    write.csv(nfd, filename)
    
    # Remove testnames rows not present in maindata
    evarcond = !(testdf$evar %in% evarnf)
    svarcond = !(testdf$svar %in% svarnf)
    mvarcond = !(testdf$mvar %in% mvarnf)
    
    testdf = testdf[evarcond & svarcond & mvarcond,]
    
    return(testdf)
  }
  
  ## Dichotomoize SNP
  dichot_df = function(df){
    df = data.frame(apply(df,2, function(x){
      y = table(x)
      if(length(y) >2){
        if(y[1] > y[3]){z = ifelse(x == 0,0,2)}
        else if(y[1] < y[3]){z = ifelse(x == 2,2,0)}
        else{z = ifelse(x == 0,0,2)}
      }
      else{z = x}
      z
    }))
    a = lapply(1:ncol(df), function(x){
      z = df[, x]
      if(any(z == 1)){
        if(any(z==0)){z= ordered(z, levels = c(0,1))}
        else{z= ordered(z, levels = c(1,2))}
      }
      else{z= ordered(z, levels = c(0,2))}
      z = data.frame(z=z)
    })
    a = do.call(cbind,a)
    colnames(a) = names(df)
    df = a
    return(df)
  }
  
  ####################
  ## Generate dataset
  ###################
  
  datagen = function(testdf, rownum, sampledf, beta = c(0,0,1)){
    y = unlist(testdf[rownum,]);
    
    # Get the real data
    df=data.frame(evar = sampledf[[1]][,y[1]], svar = sampledf[[2]][,y[2]], mvar = sampledf[[3]][,y[3]], stringsAsFactors = F)
    names(df) = c("evar", "svar", "mvar")
    
    ## generating randome values 
    set.seed(2)
    ## generating 128 random error 
    random_value = rnorm(n=nrow(df), mean=0, sd=0.25)
    
    idf = df[,2:3]    ## only SNP and methylation 
    idf$int = idf[,1]*idf[,2]    ## SNP*methylation 
    
    edf = data.frame(mapply(`*`,idf[,c("svar","mvar", "int")],beta))  ## multiple
    
    edf$con = 1
    edf$rand = random_value
    
    ## artificial gene expression y = beta0 + beta1*SNP + beta2*meCpG + beta3*SNP x meCpG + error
    edf$tot = rowSums(edf)
    
    ## replacing  real gene expression with artifical values 
    df$evar = mosaic::zscore(edf$tot)
    
    return(df)
  }
  
  ## Test
  mfit = function(df){
    # Model 1
    m1 = lm(evar ~ svar, data = df)
    # Model 2
    m2 = lm(evar ~ svar + mvar, data = df)
    # Model 3
    m3 = lm(evar ~ svar*mvar, data = df)
    
    lres13 = anova(m1, m3)
    pval13 = lres13$`Pr(>F)`[2]
    
    lres23 = anova(m2, m3)
    pval23 = lres23$`Pr(>F)`[2]
    return(list(pval13, pval23, m1, m2, m3))
  }
  
  
  modelfeat = function(mfit){
    
    # Model pvalue
    fstat = summary(mfit)$fstatistic
    pval = pf(fstat[1], fstat[2], fstat[3], lower.tail = F)
    
    # Interaction term p value
    if(nrow(summary(mfit)$coefficients) == 4){
      ipval = summary(mfit)$coefficients[4,4]
    }
    else{ipval = NA}
    
    mfeat = c(modelpval = pval, intpval = ipval)
    return(mfeat)
  }
  
}

###########################
##########################
#### Load main file ####
{
  Prefile = cohortdf(fileloc=".", pattern = "filtered.txt")
  
  #### Prepare Models 1 to 3 ####
  # Prepare the genomic models list
  testnames = genedataprep(fileurl = "./cpg_snp_gene_combination_for_memo-eQTL_mapping.txt")
  testnames = cleangenedata(maindf = Prefile, testdf = testnames, savename = "full")
  
  # Remove SNPs not present in pruning list
  filtersnp <- read.delim("./CPGEA_gdata_rsID_filtered_LD_pruned_0.8.txt", header = F)
  testnames = testnames[testnames$svar %in% filtersnp$V1,]
  
  # Normalize the SNP in range 0,1
  sampledf = dataprep(maindf = Prefile, testdf = testnames)
  originalsampledf = sampledf
  sampledf[[2]] = sampledf[[2]]/2
  sampledf[[1]] = apply(sampledf[[1]],2,RNOmni::RankNorm)
  
  # Randomly remove Identical SNPs
  snpdf = sampledf[[2]]
  
  # Prepare the pairs of identical SNPs
  snpgenedt = setDT(data.frame(t(snpdf)))
  snpduplist = snpgenedt[snpgenedt[, .N, by=names(snpgenedt)][ N > 1L ], on=names(snpgenedt), .(N, locs = .(.I)), by=.EACHI]
  int_pairs = lapply(1:nrow(snpduplist), function(x) {#
    snps = unlist(snpduplist$locs[x])
    data.frame(t(data.frame(names(snpdf)[snps])))
  })
  int_pairs = do.call(rbind.fill,int_pairs)
  
  # Find duplicates (gene and meth) in testnames
  methgenedt = setDT(testnames[,-2])
  duplist = methgenedt[methgenedt[, .N, by=names(methgenedt)][ N > 1L ], on=names(methgenedt), .(N, locs = .(.I)), by=.EACHI]
  
  foculist = duplist
  plan(multisession(workers = 1))
  
  droprowlist = lapply(1:nrow(foculist), function(x){#future.apply::future_
    r = unlist(foculist$locs[x])
    svarname = testnames$svar[r]
    tempsnp = snpdf[,svarname]
    
    rescor = tempsnp[!duplicated(as.list(tempsnp))]
    if(ncol(tempsnp) == ncol(rescor)){res = NA}
    else{
      if(ncol(rescor) == 1){res = r[-1]}
      else{res = r[-which(svarname %in% names(rescor))]}
    }
  })
  droprowlist = unlist(droprowlist)
  droprowlist = droprowlist[!is.na(droprowlist)]
  testnames = testnames[-droprowlist,]
  rownames(testnames) = 1:nrow(testnames)
  
  ### for testing 
  #testnames_ori <- testnames
  #testnames <- testnames_ori[1:300, ]    ## testing 300 lines
  
  #### Create Artificial Dataset ####
  betalist = list()
  betalist[[1]] = c(0,0,0)
  betalist[[2]] = c(0,0,0.1)
  betalist[[3]] = c(0,0,0.5)
  betalist[[4]] = c(0.1,0,0)
  betalist[[5]] = c(0.1,0.1,0.1)
  betalist[[6]] = c(0.1,0.1,-0.1)
  betalist[[7]] = c(0.5,0.5,-0.5)
  betalist[[8]] = c(0.1,0.1,0.5)
  betalist[[9]] = c(0.5,0.5,0.1)
  betalist[[10]] = c(0.1,0.1,-0.5)
  betalist[[11]] = c(0.5,0.5,-0.1)
  betalist[[12]] = c(0.5,0.5,0.5)
  #betalist[[13]] = c(0.5,0,0)
  
  plan(multisession(workers = 4))
  res = future.apply::future_lapply(1:nrow(testnames), function(x){
    if(x%%length(betalist) == 0){beta = betalist[[length(betalist)]]}
    else{beta = betalist[[x%%length(betalist)]]}
    datagen(testdf = testnames, rownum = x, sampledf = sampledf, beta = beta)
  }, future.seed = T)
  
  plan(multisession(workers = 1))
  mfit = future.apply::future_lapply(res, function(df){fit = mfit(df = df)}, future.seed = T)
  
  
  ## 
  adf = lapply(1: length(mfit), function(x){
    l = mfit[[x]]
    m13 = l[[1]]; m23 = l[[2]]
    
    if(any(is.na(c(m13,m23)))){
      print(x)
      df = data.frame(m13pval = m13, m23pval = m23, m1pval = NA, m2pval = NA, m3pval = NA, m3intpval = NA)
    }
    else{
      m1 = modelfeat(mfit = l[[3]])
      m2 = modelfeat(mfit = l[[4]])
      m3 = modelfeat(mfit = l[[5]])
      df = data.frame(m13pval = m13, m23pval = m23, m1pval = m1[1], m2pval = m2[1], m3pval = m3[1], m3intpval = m3[2])
    }
    
    df = cbind(testnames[x,], df)
    df$num = x
    
    ## if beta[3] == 0, there is no interaction  
    ## circulating 12 combinations 
    if(x%%12 == 0){beta = betalist[[12]]}
    else{beta = betalist[[x%%12]]}
    
    df$snp = beta[1]
    df$cpg = beta[2]
    df$snp_cpg = beta[3]
    df$true = ifelse(beta[3] == 0,0,1)
    df
  })
  
  adf = do.call(rbind, adf)
  
  saveRDS(adf, "Simulation_results.RDS")
  
  #### RUN FDR ANALYSIS ####
  ## check the interaction only 
  fdrdf = adf[!is.na(adf$m1pval),]
  fdrpval = p.adjust(fdrdf$m3intpval, method = "fdr", n = nrow(fdrdf))
  fdrdf$predpval = fdrpval
  fdrdf$pred = ifelse(fdrpval < 0.05, 1, 0)
  
  MLmetrics::ConfusionMatrix(y_pred = fdrdf$pred, y_true = fdrdf$true)
  MLmetrics::F1_Score(y_pred = fdrdf$pred, y_true = fdrdf$true) # F1 Score = 0.305
  ## true positive : 1165
  
  #### Proposed Method Analysis ####
  prdf = adf[!is.na(adf$m1pval),] # number of rows: 117474
  prdf$pred = ifelse(prdf$m3pval < 0.05 & prdf$m13pval < 0.05 & prdf$m23pval < 0.05, 1, 0)
  MLmetrics::ConfusionMatrix(y_pred = prdf$pred, y_true = prdf$true)
  MLmetrics::F1_Score(y_pred = prdf$pred, y_true = prdf$true) # F1 Score = 0.333
  
  ## true positive : 2871
}
