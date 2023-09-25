#####################################################################
## run memo-eQTL mapping, which will take ~24 hours for full data set

## for testing: top 150 gene-snp-cpg were tested in this script 
## for full dataset needs to set:
## res_nrow <- nrow(testnames)

############################################################################

rm(list = ls())
setwd("./Results")

###########################
#### Load the packages ####
library(readr)        ## read_table2
library(data.table)
library(plyr)
library(future.apply)
library(future)
library(dplyr)
library(distributions3)


###############################
#### Define the functions ####
##############################
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
}

## Genomic dataset preparation
genedataprep = function(fileurl){
  
  # Get the names of the columns provided by Yong
  testfile <- read_table2(fileurl, col_names = F)
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
  #write.csv(nfd, filename)
  
  svarnf = setdiff(testdf$svar, names(sampledf[[2]]))
  nfd = data.frame(snp = svarnf)
  filename = paste("snpnf",savename,".csv", sep="")
  #write.csv(nfd, filename)
  
  mvarnf = setdiff(testdf$mvar, names(sampledf[[3]]))
  nfd = data.frame(meth = mvarnf)
  filename = paste("methnf",savename,".csv", sep="")
  #write.csv(nfd, filename)
  
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

########################################
# Model building and analysis functions
## Bonferroni Test
bonferroni = function(a=0.05,k=100) return(1-((1-a)^(1/k)))

## Compare two model performance
modelperf = function(evar, svar, mvar, df =NA, dflist = sampledf, combtype = c("m1m3", "m2m3", "m1m2")){
  if(sum(is.na(df))){
    # Get the data
    df=data.frame(evar = dflist[[1]][,evar], svar = dflist[[2]][,svar], mvar = dflist[[3]][,mvar], stringsAsFactors = F)
    # str(df)
    names(df) = c("evar", "svar", "mvar")
  }
  
  if(combtype == "m1m3"){
    # Model 1
    m1 = lm(evar ~ svar, data = df)
    # Model 2
    m2 = lm(evar ~ svar*mvar, data = df)
  }
  else if(combtype == "m2m3"){
    # Model 1
    m1 = lm(evar ~ svar + mvar, data = df)
    # Model 2
    m2 = lm(evar ~ svar*mvar, data = df)
  }
  else{
    # Model 1
    m1 = lm(evar ~ svar, data = df)
    # Model 2
    m2 = lm(evar ~ svar + mvar, data = df)
  }
  lres = anova(m1, m2)
  
  pval = lres$`Pr(>F)`[2]
  return(list(pval, m1, m2))
}

## Extract model features
modelfeat = function(model){
  # print(summary(model)$coefficients)
  svar = summary(model)$coefficients[2,c(1,4)]
  mvar = summary(model)$coefficients[3,c(1,4)]
  r2 = summary(model)$adj.r.squared
  fstat = summary(model)$fstatistic
  pval = pf(fstat[1], fstat[2], fstat[3], lower.tail = F)
  modelval = c(svar = svar, mvar=mvar, rsquare = r2, pval = pval)
  return(modelval)
}


basicmodelfeat = function(model){
  modsum = summary(model)
  ## Number of Features
  numfeat = nrow(modsum$coefficients)-1
  modelconfint = confint(model)
  
  # Get beta and 95% CI of each variable
  varmetric = lapply(1:numfeat, function(x){
    featval = modsum$coefficients[x+1,]
    featname = rownames(modsum$coefficients)[x+1]
    featbeta = featval[1]
    featpval = featval[4]
    featlci = modelconfint[x+1,1]
    featuci = modelconfint[x+1,2]
    c(featname = featname, featbeta = featbeta, featpval = featpval, featlci = featlci, featuci = featuci)
  })
  
  allvarmetric = do.call(cbind, varmetric)
  
  # Get model R2 and pvalue
  modelrsq = modsum$adj.r.squared
  
  ## model pvalue: It is not directly given, so need to calculate from f statistic
  fstat = modsum$fstatistic
  modelpval = stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = F)
  
  modelmetric = c(modelrsq = modelrsq, modelpval = modelpval)
  all_metric = c(varmetric, modelmetric)
  all_metric = unlist(all_metric)
  return(all_metric)
}


## Run the M1, M2 and M3 models
three_modelperf = function(evar, svar, mvar, df = NA, dflist = sampledf){
  #if(is.na(df)){
  if(sum(is.na(df))){
    # Get the data
    df=data.frame(evar = dflist[[1]][,evar], svar = dflist[[2]][,svar], mvar = dflist[[3]][,mvar], stringsAsFactors = F)
    # str(df)
    names(df) = c("evar", "svar", "mvar")
  }
  
  # Model 1
  m1 = lm(evar ~ svar, data = df)
  
  # Model 2
  m2 = lm(evar ~ svar + mvar, data = df)
  # Model 2
  m3 = lm(evar ~ svar*mvar, data = df)
  lres13 = anova(m1, m3)
  pval13 = lres13$`Pr(>F)`[2]
  
  lres23 = anova(m2, m3)
  pval23 = lres23$`Pr(>F)`[2]
  return(list(pval13, pval23, m1, m2, m3))
}


outrem_loopperf = function(testdf, rownum, sampledf){
  y = unlist(testdf[rownum,]);
  
  # Get the data
  df=data.frame(evar = sampledf[[1]][,y[1]], svar = sampledf[[2]][,y[2]], mvar = sampledf[[3]][,y[3]], stringsAsFactors = F)
  names(df) = c("evar", "svar", "mvar")
  # print(table(df$svar))
  zval = scale(df$evar)[,1]
  # plot(density(zval))
  outrow = which(zval > 4 | zval < -4)
  # print(outrow)
  if(length(outrow)>0){df = df[-outrow, ]}
  
  if(length(table(df$svar)) <2){out=c(m1m3pval = NA, m2m3pval = NA); return(out)}
  
  modelres = three_modelperf(evar = y[1], svar = y[2], mvar = y[3], dflist = sampledf, df = df)
  
  m1m3pval = modelres[[1]]
  m2m3pval = modelres[[2]]
  # print(modelres)
  sumfit = summary(modelres[[5]])$coefficients[,c(1,4)]
  # print(summary(modelres[[3]]))
  # Get snp:meth interaction names
  intnum = grep(":", rownames(sumfit))
  sumfit = data.frame(sumfit[intnum,])
  
  snp1beta = sumfit[1,1]; snp1pval = sumfit[2,1]
  intres = c(snp1beta = snp1beta, snp1pval = snp1pval)
  
  modelval2 = modelfeat(model = modelres[[4]])
  modelval3 = modelfeat(model = modelres[[5]])
  
  model1metrics = basicmodelfeat(modelres[[3]])
  model2metrics = basicmodelfeat(modelres[[4]])
  model3metrics = basicmodelfeat(modelres[[5]])
  out = c(m1m3pval = m1m3pval,m2m3pval = m2m3pval, intres, modelval2, modelval3, m1m = model1metrics, m2m = model2metrics, m3m = model3metrics)
  return(out)
}

## Regression Clustering
rc_clusdataprep = function(df, outlier_rem = T, Dichot = T){
  if(outlier_rem){
    zval = scale(df$evar)[,1]
    # plot(density(zval))
    outrow = which(zval > 4 | zval < -4)
    # print(outrow)
    if(length(outrow)>0){df = df[-outrow, ]}
  }
  if(Dichot){
    y = table(df$svar)
    if(length(y)>2){
      if(y[1] >= y[3]){
        if(min(y) == y[3]){df$svar[df$svar == 2] = 1}
        else{df$svar[df$svar == 1] = 2}
      }
      else{
        if(min(y) == y[1]){df$svar[df$svar == 0] = 1}
        else{df$svar[df$svar == 1] = 0}
      }
    }
  }
  newdf = df
  newdf$svar = ordered(newdf$svar, levels = c(min(newdf$svar), max(newdf$svar)))
  return(list(df = df, newdf = newdf))
}

svarmodelpara = function(model){
  coef = coef(model)[2]
  ci = confint(model)[2,]
  pval = summary(model)$coefficients[2,4]
  rsquare = summary(model)$adj.r.squared
  fstat = summary(model)$fstatistic
  modelpval = pf(fstat[1], fstat[2], fstat[3], lower.tail = F)
  modelpara = c(coef = coef, ci = ci, pval = pval, modelpval = modelpval, rsquare = rsquare)
  return(modelpara)
}

regclusres = function(df, k=2, cutpoint = NA, rownumber = 1:10){
  if(is.na(cutpoint)){
    df1 = df[rownumber,]
    df2 = df[-rownumber,]
  }
  else{
    df1 = df[df$mvar < cutpoint,]
    df2 = df[df$mvar >= cutpoint,]
  }
  # print(table(df1$svar))
  # print(table(df2$svar))
  
  if(nrow(df1) < 0.25*nrow(df) | nrow(df2) < 0.25*nrow(df)){perf = -1; svar1=svar2=res = NA; pval = 1}
  else if(any(table(df1$svar) == 0) | any(table(df2$svar) == 0) ){ perf = -1; svar1=svar2=res = NA; pval = 1}
  else{
    l1 = lm(evar~svar, df1)
    l2 = lm(evar~svar, df2)
    
    l1sum = summary(l1)$coefficient
    l2sum = summary(l2)$coefficient
    if(nrow(l1sum) <2 | nrow(l2sum) <2){perf = -1; svar1=svar2=res = NA; pval = 1}
    else{
      # print(nrow(summary(l1)$coefficient))
      # Get coefficient, se and pval
      l1para = l1sum[2,c(1,2,4)]
      l2para = l2sum[2,c(1,2,4)]
      
      z_stat = abs(l1para[1] - l2para[1])/sqrt(l1para[2]^2 + l2para[2]^2)
      Z = Normal(0,1)
      pval = 1 - cdf(Z, z_stat) + cdf(Z, -z_stat)
      perf = -log(pval); res = list(l1,l2)
      
      # Extract model data
      svar1 = svarmodelpara(model = res[[1]])
      svar2 = svarmodelpara(model = res[[2]])
    }
  }
  
  return(list(perf=perf, res= res, df1 = df1, df2=df2, pval = pval, svar1=svar1, svar2 = svar2))
}


RC = function(df,mvar = T){
  gadf = df
  gadf$evar = gadf$evar/max(gadf$evar)
  
  if(mvar){
    ga_fun = function(x){
      
      res = regclusres(df=gadf, k=2, cutpoint = x, rownumber = 1:10)
      perf = res$perf
      return(perf)
    }
    GA = GA::ga(type = "real-valued", fitness = ga_fun, lower = min(gadf$mvar), upper = max(gadf$mvar), maxiter = 100, run=10, pcrossover = 0.8, pmutation = 0.2, popSize = 100, monitor = F, parallel = F, seed = 3)
  }
  else{
    ga_fun = function(x){
      rownumber = which(x==1)
      res = regclusres(df=gadf, k=2, cutpoint = NA, rownumber = rownumber)
      perf = res$perf
      return(perf)
    }
    GA = GA::ga(type = "binary", fitness = ga_fun, nBits = nrow(gadf), maxiter = 100, run=10, pcrossover = 0.8, pmutation = 0.2, popSize = 100, monitor = F, parallel = F, seed = 3)
  }
  return(GA)
}


svarlabel = function(svarpara){
  betalab = paste(round(svarpara[1],2), "(",round(svarpara[2],2), ",", round(svarpara[3],2), ")",
                  " p <",ifelse(svarpara[4] < 0.001, "0.001", round(svarpara[4],3)), sep = "")#,,
  modellab = paste("model = p<", ifelse(svarpara[5] < 0.001, "0.001", round(svarpara[5],3)), " rsq = ", round(svarpara[6],2), sep = "")
  return(list(betalab, modellab))
}


## model feature property
modelprop_loopperf = function(testdf, rownum, sampledf){
  y = unlist(testdf[rownum,]);
  
  # Get the data
  df=data.frame(evar = sampledf[[1]][,y[1]], svar = sampledf[[2]][,y[2]], mvar = sampledf[[3]][,y[3]], stringsAsFactors = F)
  names(df) = c("evar", "svar", "mvar")
  # print(table(df$svar))
  zval = scale(df$evar)[,1]
  # plot(density(zval))
  outrow = which(zval > 4 | zval < -4)
  # print(outrow)
  if(length(outrow)>0){df = df[-outrow, ]}
  
  # Sample size
  ss = nrow(df)
  # mean and sd
  gene_mean = mean(df$evar)
  gene_sd = sd(df$evar)
  meth_mean = mean(df$mvar)
  meth_sd = sd(df$mvar)
  # svar distribution
  svartable = table(df$svar)
  if(any(names(svartable) %in% c(0))){svar_0 = svartable[names(svartable) == c(0)]}
  else{svar_0 = 0}
  
  if(any(names(svartable) %in% c(0.5))){svar_1 = svartable[names(svartable) == c(0.5)]} 
  else{svar_1 = 1}
  
  if(any(names(svartable) %in% c(1))){svar_2 = svartable[names(svartable) == c(1)]} 
  else{svar_2 = 0}
  
  svar_ss = sum(svartable)
  
  out = c(ss = ss, gene_mean = gene_mean, gene_sd = gene_sd, meth_mean = meth_mean, meth_sd = meth_sd,
          svar_0 = svar_0, svar_1 = svar_1, svar_2 = svar_2, svar_ss = svar_ss)
  return(out)
}
# Randomly remove Identical SNPs

fastCorrel=function(inputdf, cutoff){
  f = function(df, cut = cutoff){
    dfcor=cor(df, use = "pairwise.complete.obs", method = "spearman")
    highcor=apply(dfcor,1,function(x) {y = length(which(x>=cut | x<= (-1*cut))); y})
    
    var=attributes(highcor)$names
    
    # Select variables with no correlation
    selvar = names(df)[which(highcor <= 1)]
    
    # Keep variables with high correlation
    highvar = setdiff(var, selvar)
    
    varseq = data.frame(variable = attributes(highcor)$names, ncor = highcor)
    varseq = varseq[varseq$ncor >1, ]
    varseq = varseq[order(varseq$ncor, decreasing = T),]
    
    # print(list(selvar,varseq))
    return(list(selvar,varseq))
  }
  
  outf = f(df = inputdf, cut = cutoff)
  selvar = outf[[1]]
  varseq = outf[[2]]
  
  for(i in 1:nrow(varseq)){
    if(nrow(varseq) <3){
      selvar = c(selvar, varseq$variable[-1])
      break}
    newvar = varseq$variable[-1]
    newdf = inputdf[,newvar]
    outf = f(df = newdf)
    selvar = c(selvar,outf[[1]])
    varseq = outf[[2]]
  }
  
  finaldf = inputdf[,selvar]
  return(finaldf)
}
}

#######################################################
#### Load main files and perform memo-eQTL mapping ####
#######################################################
{
############################
### load data and filtering 
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
sampledf[[1]] = apply(sampledf[[1]],2,RNOmni::RankNorm)   ## apply RankNorm for gene expression

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
dim(testnames)  
## 90959     3   ## after further filtering 
}

#########################################
# Prepare the models for each combination
{
plan(multisession(workers = 1))

## reuslts for the all combinations !!!!!!
# res_nrow <- nrow(testnames)
## testing first 100 records
res_nrow = 150
res = future.apply::future_lapply(1:res_nrow, function(x){
  outrem_loopperf(testdf = testnames, rownum = x, sampledf = sampledf)})

# Combine the results
reslist = do.call(bind_rows, res)
normalres = data.frame(testnames[1:res_nrow, 1:3], data.frame(reslist[ ,c(1:2,17:52)]))

varnames = c("m1_svar_beta", "m1_svar_pval", "m1_svar_lci", "m1_svar_uci","m1_model_rsq","m1_model_pval",
             "m2_svar_beta", "m2_svar_pval", "m2_svar_lci", "m2_svar_uci",
             "m2_mvar_beta", "m2_mvar_pval", "m2_mvar_lci", "m2_mvar_uci",
             "m2_model_rsq","m2_model_pval",
             "m3_svar_beta", "m3_svar_pval", "m3_svar_lci", "m3_svar_uci",
             "m3_mvar_beta", "m3_mvar_pval", "m3_mvar_lci", "m3_mvar_uci",
             "m3_svar_mvar_beta", "m3_svar_mvar_pval", "m3_svar_mvar_lci", "m3_svar_mvar_uci",
             "m3_model_rsq","m3_model_pval")

removefeaturenames = grepl("featname", names(normalres))
normalres = normalres[,!removefeaturenames]
names(normalres)[-1:-5] =  varnames
m1to3 = normalres
}

###################################
#### Perform cluster Analysis ####
# Prepare the genomic models list
{
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

foculist = duplist#[duplist$evar == "CD320" & duplist$mvar == "chr19:8347912",]
plan(multisession(workers = 1))

droprowlist = lapply(1:nrow(foculist), function(x){#future.apply::future_
  # cat(x, " ")
  r = unlist(foculist$locs[x])
  svarname = testnames$svar[r]
  tempsnp = snpdf[,svarname]
  
  #rescor = fastCorrel(inputdf = tempsnp, cutoff = 1)
  # rescor = data.frame(rescor)
  rescor = tempsnp[!duplicated(as.list(tempsnp))]
  # str(rescor)
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

## Get Sampledf
binsampledf = dataprep(maindf = Prefile, testdf = testnames)
binsampledf[[1]] = apply(binsampledf[[1]],2,RNOmni::RankNorm)

fulldflist = list()
nrow(testnames)

######################################
# Perform M1M2 and Optimal clustering
testnames <- testnames[1:res_nrow, ]    ## for testing  

plan(multisession(workers = 8))
cutpoint = 100           ## each block has 100 tests
nsmart = floor(nrow(testnames)/cutpoint)
smartfull = pbapply::pblapply(1:(nsmart+1), function(smart){
  cat("Block ", smart)
  lb = 1+ cutpoint*(smart-1)
  ub = min(cutpoint + cutpoint*(smart-1), nrow(testnames))
  
  fullres = future.apply::future_lapply(lb:ub, function(i){
    
    # cat(i, " ")
    y = unlist(testnames[i,]);
    
    # Get the data
    df=data.frame(evar = sampledf[[1]][,y[1]], svar = sampledf[[2]][,y[2]], mvar = sampledf[[3]][,y[3]], stringsAsFactors = F)
    # str(df)
    
    # m1m2
    res = modelperf(df =df, dflist = sampledf, combtype = c("m1m2"))
    m1m2pval = res[[1]]
    # print(m1m2pval)
    
    # Do clustering
    mod_df=data.frame(evar = binsampledf[[1]][,y[1]], svar = binsampledf[[2]][,y[2]], mvar = binsampledf[[3]][,y[3]], stringsAsFactors = F)
    names(mod_df) = c("evar", "svar", "mvar")
    if(length(unique(mod_df$svar)) == 1){
      outdata = c(y, m1m2pval = m1m2pval)
    }
    else{
      dflist = rc_clusdataprep(df = mod_df, outlier_rem = T, Dichot = T)
      # str(dflist)
      mod_df = dflist$newdf
      # print(T)
      # Regression Clustering
      res = RC(df=mod_df)
      # str(res)
      finres = regclusres(df = mod_df, cutpoint = mean(res@solution))
      
      df1 = finres$df1
      clus1 = paste("Low",round(min(df1$mvar),2),round(mean(df1$mvar),2),round(max(df1$mvar),2),sep = "_")
      df1$cluster = clus1
      df2 = finres$df2
      clus2 = paste("High",round(min(df2$mvar),2),round(mean(df2$mvar),2),round(max(df2$mvar),2), sep = "_")
      df2$cluster = clus2
      
      outdata = c(y, m1m2pval = m1m2pval, clus1label = clus1, clus1 = finres$svar1, clus2label = clus2, clus2 = finres$svar2)
      
    }
    outdata = data.frame(outdata, stringsAsFactors = F)
    outdata = data.frame(t(outdata), stringsAsFactors = F)
    
  }, future.seed = T)
  
  out = do.call(rbind.fill, fullres)
  # filename = paste(smart, "m1m2_cluster.csv", sep = "")
  # write.csv(out, filename)
  out
})
csvdf = lapply(smartfull, function(x) if(ncol(x) ==18){x}else{x[,1:18]})
clusterdf = do.call(rbind, csvdf)
}

###########################
#### Prepare model 4 ####
{
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
sampledf[[3]] = data.frame(apply(sampledf[[3]],2,function(x)scale(x)[,1]))
names(sampledf[[3]]) = stringi::stri_replace_all_fixed(names(sampledf[[3]]), ".", ":")

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

foculist = duplist#[duplist$evar == "CD320" & duplist$mvar == "chr19:8347912",]
plan(multisession(workers = 3))

droprowlist = lapply(1:nrow(foculist), function(x){#future.apply::future_
  # cat(x, " ")
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

## Get Sampledf
binsampledf = dataprep(maindf = Prefile, testdf = testnames)
binsampledf[[1]] = apply(binsampledf[[1]],2,RNOmni::RankNorm)

fulldflist = list()
nrow(testnames)

## for M4
testnames <- testnames[1:res_nrow, ]    ## for testing   

cutpoint = 100
nsmart = floor(nrow(testnames)/cutpoint)
smartfull = pbapply::pblapply(1:(nsmart+1), function(smart){#
  cat("Block ", smart)
  lb = 1+ cutpoint*(smart-1)
  ub = min(cutpoint + cutpoint*(smart-1), nrow(testnames))
  fullres = lapply(lb:ub, function(i){#future.apply::future_ 
    
    y = unlist(testnames[i,]);
    
    # Get the data
    df=data.frame(evar = sampledf[[1]][,y[1]], svar = sampledf[[2]][,y[2]], mvar = sampledf[[3]][,y[3]], stringsAsFactors = F)
    
    if(length(unique(df$svar)) == 1){
      outdata = data.frame(svar_percent =NA, mvar_percent = NA, svar.mvar_percent = NA, m3r2 = NA, m4svarbeta =NA, m4svarpval = NA, m4r2 = NA, stringsAsFactors = F)
    }
    else{
      # get individual model contribution in moderation model
      fit1 = lm(evar ~ svar*mvar, data = df)
      if(any(is.na(fit1$coefficients))){
        outdata = data.frame(svar_percent =NA, mvar_percent = NA, svar.mvar_percent = NA, m3r2 = NA, m4svarbeta =NA, m4svarpval = NA, m4r2 = NA, stringsAsFactors = F)
      }
      else{
        # print(fit1)
        metrics <-relaimpo::calc.relimp(fit1, type = c("lmg"), rela = T)
        fit1r2 = summary(fit1)$adj.r.squared
        # print(metrics)
        fit2 = lm(mvar ~ svar, data = df)
        linmodsum = summary(fit2)
        beta_pval = linmodsum$coefficients[2,c(1,4)]
        fit2var = c(m4svarbeta = beta_pval[1], m4svarpval = beta_pval[2], r2 = linmodsum$r.squared)
        # print(summary(fit2))
        
        outdata = data.frame(metrics@lmg, stringsAsFactors = F)
        outdata = data.frame(t(outdata), stringsAsFactors = F)
        names(outdata) = paste(names(outdata),"percent", sep = "_")
        outdata[,c("m3r2", "m4svarbeta","m4svarpval","m4r2")] = c(fit1r2, fit2var)
      }
    }
    outdata[,c("gene", "snp", "cpg")] = testnames[i,]
    outdata
  })#, future.seed = T
  out = do.call(rbind.fill, fullres)
  out = data.frame(out)
  out
})

csvfiles = lapply(smartfull, function(x) if(ncol(x) ==10){x}else{x[,1:10]})
csvdf = do.call(rbind, csvfiles)
m4 = csvdf
}

###########################################
#### Pool data from all three analysis ####
df_fin = merge(m1to3, clusterdf, by.x = c('evar', 'svar', 'mvar'), by.y = c('evar', 'svar', 'mvar'))
df_fin = merge(df_fin, m4, by.x = c('evar', 'svar', 'mvar'), by.y = c('gene', 'snp', 'cpg'))

#saveRDS(df_fin, "all_tested_snp_cpg_gene_cmb_4_memo-eqtl.RDS")
saveRDS(df_fin, "150_tested_snp_cpg_gene_cmb_4_memo-eqtl.RDS")

}
