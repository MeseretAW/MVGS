#library
library(readxl)
library(gaston)
library(tidyverse)
library(asreml)
library(psych)
library(reshape2)
library(ggcorrplot)
library(ggplot2)

#set working directory
setwd("MVGS_for_quality_traits")


## ------------------------------------------------------------------------------------------------------------
#read in phenotypic data
pheno <- read.table("Phenotypic_data.txt", header = TRUE, sep = "\t") # Change sep to your delimiter

pheno<-pheno %>%
  mutate(env=paste(YR, EXPT, NURNAME, sep = "_")) %>%
  select(env,everything())

pheno[,1:10]<-lapply(pheno[,1:10], as.factor)
pheno[,11:ncol(pheno)]<-lapply(pheno[,11:ncol(pheno)], as.numeric)

#read genotypic data and calculate GRM
geno <- read.table("Genotypic_data.txt", header = TRUE, sep = "\t") # Change sep to your delimiter

geno_amat<-GRM(geno, autosome.only = F)


## ------------------------------------------------------------------------------------------------------------
#making an object for results
blues<-pheno %>%
  distinct(ID)

#define traits
traits<-colnames(pheno)[11:ncol(pheno)]

#loop
for(i in traits){
  
  #print
  print(paste("Performing Multi-Environmental Analysis for", i))
  
  #pull out the data
  a<-pheno %>% 
    select(ID, ROW, COLUMN, env, all_of(i)) %>%
    drop_na()
  colnames(a)[ncol(a)]="y"
  a[,1:4]<-lapply(a[,1:4], as.factor)
  a$y<-as.numeric(a$y)
  
  #run the model
  fit<-asreml(fixed = y~ID,
              random = ~ idv(ROW):idv(env)+idv(COLUMN):idv(env)+idv(env),
              residual = ~ idv(units),
              data = a,
              trace = F)
  
  #pull estimates
  pred<-predict(fit, classify = "ID")[["pvals"]][,1:3]
  
  #rename columns
  colnames(pred)[2:3]=c(i, paste("SE_", i, sep = ""))
  
  #join
  blues<-blues %>%
    left_join(pred, by="ID")
  
  #remove
  remove(a, fit, pred)
}

#select traits of interest
blues_of_interest<-blues %>% 
  select(-colnames(blues)[grep("SE_", colnames(blues))])%>% 
  select(ID,
         SRCW,
         Hard,
         Dia,
         QuadT,
         QuadBC)

#produce correlation matrix
png("rplot.png", width = 10, height = 10, units = "in", res = 300)

pairs.panels(blues_of_interest[,-1], 
             stars = T,
             lm = T,
             hist.col = "gray")

dev.off()

#select highly correlated traits
blues_cor<-cor(blues_of_interest[,-1], use = "complete.obs")

#linearize the corr matrix
melt_cor<-melt(blues_cor)

#selecting highly correlated traits
SRCW_cor<-melt_cor %>% 
  filter(Var2=="SRCW") %>% 
  filter(value>=0.25 |
           value<=-0.25)



## ------------------------------------------------------------------------------------------------------------
#make GBS data into a matrix
geno<-as.matrix(geno)

#do the pca
pca<-prcomp(geno)
summary(pca)$importance[,1:5]
pca<-as.data.frame(pca$x)

#plotting the pca
plot(pca$PC1,pca$PC2)
plot(pca$PC2,pca$PC3)
plot(pca$PC1, pca$PC3)


## ------------------------------------------------------------------------------------------------------------
#format the data
blues_of_interest$ID<-as.factor(blues_of_interest$ID)
blues_of_interest[,2:ncol(blues_of_interest)]<-lapply(blues_of_interest[,2:ncol(blues_of_interest)],as.numeric)

#subset individuals in the GRM found in my data
geno_amat<-geno_amat[rownames(geno_amat) %in% blues_of_interest$ID,
                     colnames(geno_amat) %in% blues_of_interest$ID]

#select only genotypes found in the GRM
blues_of_interest<-blues_of_interest %>%
  filter(ID %in% rownames(geno_amat))

#write a multivariate variate model to estimate genetic correlation
fit<-asreml(fixed = cbind(SRCW,
                          Hard,
                          Dia,
                          QuadT,
                          QuadBC)~trait,
            random = ~vm(ID, geno_amat):us(trait),
            residual = ~id(units):us(trait),
            data = blues_of_interest,
            na.action = na.method(y='include', x='include'),
            workspace = "8gb",
            maxit = 50,
            trace = TRUE)

vcov<-summary(fit)$varcomp %>%
  rownames_to_column(var = "comp") %>%
  rownames_to_column(var = "xform") %>%
  mutate(xform=paste("V",xform,sep = ""))

h2_SRCW<-vpredict(fit,h2~V1/(V1+V17))
h2_Hard<-vpredict(fit,h2~V3/(V3+V19))
h2_Dia<-vpredict(fit,h2~V6/(V6+V22))
h2_QuadT<-vpredict(fit,h2~V10/(V10+V26))
h2_QuadBC<-vpredict(fit,h2~V15/(V15+V31))

h2<-rbind(h2_SRCW, h2_Hard, h2_Dia, h2_QuadT, h2_QuadBC)

mv_h2<-data.frame(Trait=c("SRCW","Hard","Dia","QuadT","QuadBC"),
                  h2=h2$Estimate,
                  SE=h2$SE)
summary(fit)
#multivarite genetic correlations
SRCW_Hard<-vpredict(fit, ~ V2 / sqrt(V1 * V3))$Estimate
SRCW_Dia<-vpredict(fit, ~ V4 / sqrt(V1 * V6))$Estimate
SRCW_QuadT<-vpredict(fit, ~ V7 / sqrt(V1 * V10))$Estimate
SRCW_QuadBC<-vpredict(fit, ~ V11 / sqrt(V1 * V15))$Estimate
Hard_Dia<-vpredict(fit, ~ V5 / sqrt(V3 * V6))$Estimate
Hard_QuadT<-vpredict(fit, ~ V8 / sqrt(V3 * V10))$Estimate
Hard_QuadBC<-vpredict(fit, ~ V12 / sqrt(V3 * V15))$Estimate
Dia_QuadT<-vpredict(fit, ~ V9 / sqrt(V6 * V10))$Estimate
Dia_QuadBC<-vpredict(fit, ~ V13 / sqrt(V6 * V15))$Estimate
QuadT_QuadBC<-vpredict(fit, ~ V14 / sqrt(V10 * V15))$Estimate


cormat<-matrix(data=c(1,SRCW_Hard,SRCW_Dia,SRCW_QuadT,SRCW_QuadBC,
                      0, 1, Hard_Dia,Hard_QuadT,Hard_QuadBC,
                      0, 0, 1, Dia_QuadT,Dia_QuadBC,
                      0, 0 , 0, 1, QuadT_QuadBC,
                      0, 0, 0, 0, 1),
               
               ncol = 5,
               nrow = 5,
               byrow = T)
cormat[lower.tri(cormat)]=cormat[upper.tri(cormat)]
colnames(cormat)<-c("SRCW","Hard","Dia",	"QuadT",	"QuadBC")
rownames(cormat)<-colnames(cormat)

#write a loop for univariate heritability
traits<-colnames(blues_of_interest)[-1]

h2<-c()

for(i in traits){
  
  #pull data
  a<-blues_of_interest[,c("ID",i)]
  a<-a[!is.na(a[,2]),]
  colnames(a)[2]="y"
  
  #pull only lines found in a
  b<-geno_amat[rownames(geno_amat) %in% a$ID,
               colnames(geno_amat) %in% a$ID]
  
  #run model
  fit<-asreml(fixed = y ~ 1,
              random = ~ vm(ID, b),
              residual = ~ idv(units),
              data = a)
  
  #pull heritability
  h<-vpredict(fit, ~V1/(V1+V2))
  rownames(h)=NULL
  h$Trait=i
  h<-h %>%
    select(Trait,
           Estimate,
           SE) %>%
    rename(h2=Estimate)
  
  #bind
  h2<-rbind(h2,h)
  
  remove(a,b,fit,h)
}


## ------------------------------------------------------------------------------------------------------------
univariate_five_fold_cv<-function(data,
                                  genotype_col_name,
                                  trait, 
                                  grm,
                                  BLUP_fixed,
                                  BLUP_random,
                                  BLUP_residual,
                                  workspace,
                                  graph_on,
                                  trace_on,
                                  n_maxit){
  
  #pull relevant data
  a<-data.frame(data) %>%
    select(all_of(genotype_col_name),
           all_of(trait))
  
  colnames(a)<-c("genotype_col_name",
                 "y_1")
  
  a<-a %>%
    drop_na(y_1)
  
  #make first two columns into levels
  a[,1]<-as.factor(a[,1])
  
  #make everything else numeric
  a[,2]<-as.numeric(a[,2])
  
  #select out only the lines in dataset from the grm
  grm<-grm[rownames(grm) %in% unique(a[,1]),
           colnames(grm) %in% unique(a[,1])]
  
  #make sure the lines in the dataframe are in the grm
  a<-a[a[,1] %in% rownames(grm),]
  
  #check if grm and data are same size
  if(isTRUE(unique(unique(a[,1]) %in% rownames(grm)))){
    #for testing
    #print(paste("The dataframe and grm have the same number of lines;",
    #            "proceeding with analysis for",
    #            traits[1]))
  }else{
    print(paste("Not all lines in the dataset are found in the grm;",
                "aborting analysis for",
                traits[1]))
    return(NA)
  }
  
  
  #move on after truth check
  #make a partition of the dataset randomly into training and test
  train=sample(rownames(grm),
               size = round(length(unique(a[,1]))*.8, digits = 0))
  test=rownames(grm)[!rownames(grm) %in% train] 
  
  if(isFALSE(unique(test %in% train))){
    #for testing
    #print(paste("Subdivsion of training and test successful;",
    #            "proceeding with analysis for",
    #            traits[1]))
  }else{
    print(paste("Subdivsion of training and test unsuccessful;",
                "check grm for duplicated name..."))
  }
  
  #you have to assign things to the global environment because asreml is stupid
  assign("a", a, envir = globalenv() )
  assign("grm", grm, envir = globalenv() )
  assign("BLUP_fixed", BLUP_fixed, envir = globalenv() )
  assign("BLUP_random", BLUP_random, envir = globalenv() )
  assign("BLUP_residual", BLUP_residual, envir = globalenv() )
  assign("workspace", workspace, envir = globalenv() )
  
  #calculate a dataset of adjusted means for the
  #the training data
  test_preds_actual<-a %>% 
    filter(genotype_col_name %in% test) %>%
    select(genotype_col_name, y_1) %>%
    rename(TBV=y_1)
  
  #move on after truth check
  #partition dataset
  a[a$genotype_col_name %in% test, "y_1"]=NA 
  assign("a", a, envir = globalenv() )
  
  #run ASREML gBLUP
  train_preds_predicted<-asreml(fixed = BLUP_fixed,
                                random = BLUP_random,
                                residual = BLUP_residual,
                                data = a,
                                trace = trace_on,
                                workspace = workspace,
                                maxit = n_maxit)
  
  #pull predictions for test
  train_preds_predicted<-predict(train_preds_predicted, classify = "genotype_col_name")$pvals
  train_preds_predicted<-train_preds_predicted %>%
    filter(genotype_col_name %in% test_preds_actual$genotype_col_name)  %>%
    select(genotype_col_name, predicted.value) %>%
    rename(GEBV=predicted.value)
  
  #put together results
  corr<-full_join(test_preds_actual, train_preds_predicted, by=c("genotype_col_name"))
  
  if(isTRUE(graph_on)){
    plot(corr$TBV, 
         corr$GEBV,
         xlab = "True Breeding Value",
         ylab = "BLUP",
         main = paste("GEBV vs BLUE for", trait))
  }else{
    #do nothing
  }
  
  #calculate correlation
  corr<-cor(corr[,2:ncol(corr)])[1,2]
  
  #return the correlation
  return(corr)
}

multivariate_five_fold_cv<-function(data,
                                    genotype_col_name, 
                                    traits,
                                    grm,
                                    BLUP_fixed,
                                    BLUP_random,
                                    BLUP_residual,
                                    workspace,
                                    graph_on,
                                    trace_on,
                                    n_maxit){
  
  #pull relevant data
  a<-as.data.frame(data) %>%
    dplyr::select(all_of(genotype_col_name),
                  all_of(traits))
  
  colnames(a)<-c("genotype_col_name",
                 paste("y_", 
                       seq(1,ncol(a[,2:ncol(a)])), 
                       sep = ""))
  
  a<-a %>%
    drop_na(y_1)
  
  #make first two columns into levels
  a[,1]<-as.factor(a[,1])
  
  #make everything else numeric
  a[,2:ncol(a)]<-lapply(a[,2:ncol(a)], as.numeric)
  
  #select out only the lines in dataset from the grm
  grm<-grm[rownames(grm) %in% unique(a[,1]),
           colnames(grm) %in% unique(a[,1])]
  
  #make sure the lines in the dataframe are in the grm
  a<-a[a[,1] %in% rownames(grm),]
  
  #check if grm and data are same size
  if(isTRUE(unique(unique(a[,1]) %in% rownames(grm)))){
    #for testing
    #print(paste("The dataframe and grm have the same number of lines;",
    #            "proceeding with analysis for",
    #            traits[1]))
  }else{
    print(paste("Not all lines in the dataset are found in the grm;",
                "aborting analysis for",
                traits[1]))
    return(NA)
  }
  
  #move on after truth check
  #make a partition of the dataset randomly into training and test
  train=sample(rownames(grm),
               size = round(length(unique(a[,1]))*.8, digits = 0))
  test=rownames(grm)[!rownames(grm) %in% train] 
  
  if(isFALSE(unique(test %in% train))){
    #for testing
    #print(paste("Subdivsion of training and test successful;",
    #            "proceeding with analysis for",
    #            traits[1]))
  }else{
    print(paste("Subdivsion of training and test unsuccessful;",
                "check grm for duplicated name..."))
  }
  
  #you have to assign things to the global environment because asreml is stupid
  assign("a", a, envir = globalenv() )
  assign("grm", grm, envir = globalenv() )
  assign("BLUP_fixed", BLUP_fixed, envir = globalenv() )
  assign("BLUP_random", BLUP_random, envir = globalenv() )
  assign("BLUP_residual", BLUP_residual, envir = globalenv() )
  assign("workspace", workspace, envir = globalenv() )
  
  #calculate a dataset of adjusted means for the
  #the training data
  test_preds_actual<-a %>% 
    filter(genotype_col_name %in% test) %>%
    select(genotype_col_name, y_1) %>%
    rename(TBV=y_1)
  
  #move on after truth check
  #partition dataset
  a[a$genotype_col_name %in% test, "y_1"]=NA 
  assign("a", a, envir = globalenv() )
  
  #run ASREML gBLUP
  train_preds_predicted<-asreml(fixed = BLUP_fixed,
                                random = BLUP_random,
                                residual = BLUP_residual,
                                data = a,
                                trace = trace_on,
                                workspace = workspace,
                                maxit = n_maxit)
  
  #pull predictions for test
  u<-as.numeric(train_preds_predicted$coefficients$fixed[1])
  blups<-train_preds_predicted$coefficients$random
  blups<-as.data.frame(blups) %>% 
    rownames_to_column(var="id") %>%
    separate(col = id, into = c("line", "trait"), sep = ":", remove = T) %>%
    separate(col = line, into=c("junk", "junk1", "junk2", "line"), sep = "_", remove = T) %>%
    select(line, trait, bu) %>%
    filter(trait=="trait_y_1") %>%
    select(line, bu) %>%
    rename(genotype_col_name=line,
           predicted.value=bu) %>%
    select(genotype_col_name, predicted.value) %>%
    rename(GEBV=predicted.value) %>%
    mutate(GEBV=GEBV+u) %>%
    filter(genotype_col_name %in% test_preds_actual$genotype_col_name)
  
  #put together results
  corr<-full_join(test_preds_actual, blups, by=c("genotype_col_name"))
  
  if(isTRUE(graph_on)){
    plot(corr$TBV, 
         corr$GEBV,
         xlab = "True Breeding Value",
         ylab = "GEBV",
         main = paste("GEBV vs TBV for", traits[1]))
  }else{
    #do nothing
  }
  
  #calculate correlation
  corr<-cor(corr[,2:ncol(corr)])[1,2]
  
  #return the correlation
  return(corr)
}


## ------------------------------------------------------------------------------------------------------------
#make object for the results
cv_results<-c()

#perms
perms<-100

#loop
for (i in 1:perms) {
  
  #announcement
  print(paste("Commencing Permutation =", i))
  print(Sys.time())
  
  #run univariate for SRCW
  one<-univariate_five_fold_cv(data = blues_of_interest,
                               genotype_col_name = "ID",
                               trait = c("SRCW"),
                               grm = geno_amat,
                               BLUP_fixed = formula(y_1 ~ 1),
                               BLUP_random = formula(~ vm(genotype_col_name, grm)),
                               BLUP_residual = formula(~idv(units)),
                               workspace = NULL,
                               graph_on = T,
                               trace_on = T,
                               n_maxit = NULL)
  
  two<-multivariate_five_fold_cv(data = blues_of_interest,
                                 genotype_col_name = "ID",
                                 traits = c("SRCW","Hard"),
                                 grm = geno_amat,
                                 BLUP_fixed = formula(cbind(y_1, y_2) ~ trait),
                                 BLUP_random = formula(~ vm(genotype_col_name, grm):us(trait)),
                                 BLUP_residual = formula(~id(units):us(trait)),
                                 workspace = "8gb",
                                 graph_on = T,
                                 trace_on = T,
                                 n_maxit = 75)
  
  three<-multivariate_five_fold_cv(data = blues_of_interest,
                                   genotype_col_name = "ID",
                                   traits = c("SRCW","Dia"),
                                   grm = geno_amat,
                                   BLUP_fixed = formula(cbind(y_1, y_2) ~ trait),
                                   BLUP_random = formula(~ vm(genotype_col_name, grm):us(trait)),
                                   BLUP_residual = formula(~id(units):us(trait)),
                                   workspace = "8gb",
                                   graph_on = T,
                                   trace_on = T,
                                   n_maxit = 75)
  
  four<-multivariate_five_fold_cv(data = blues_of_interest,
                                  genotype_col_name = "ID",
                                  traits = c("SRCW","QuadT"),
                                  grm = geno_amat,
                                  BLUP_fixed = formula(cbind(y_1, y_2) ~ trait),
                                  BLUP_random = formula(~ vm(genotype_col_name, grm):us(trait)),
                                  BLUP_residual = formula(~id(units):us(trait)),
                                  workspace = "8gb",
                                  graph_on = T,
                                  trace_on = T,
                                  n_maxit = 75)
  
  five<-multivariate_five_fold_cv(data = blues_of_interest,
                                  genotype_col_name = "ID",
                                  traits = c("SRCW","QuadBC"),
                                  grm = geno_amat,
                                  BLUP_fixed = formula(cbind(y_1, y_2) ~ trait),
                                  BLUP_random = formula(~ vm(genotype_col_name, grm):us(trait)),
                                  BLUP_residual = formula(~id(units):us(trait)),
                                  workspace = "8gb",
                                  graph_on = T,
                                  trace_on = T,
                                  n_maxit = 75)
  
  #run multivariate for SRCW
  six<-multivariate_five_fold_cv(data = blues_of_interest,
                                 genotype_col_name = "ID",
                                 traits = c("SRCW", "Hard", "Dia", "QuadT", "QuadBC"),
                                 grm = geno_amat,
                                 BLUP_fixed = formula(cbind(y_1, y_2, y_3, y_4, y_5) ~ trait),
                                 BLUP_random = formula(~ vm(genotype_col_name, grm):us(trait)),
                                 BLUP_residual = formula(~id(units):us(trait)),
                                 workspace = "8gb",
                                 graph_on = T,
                                 trace_on = T,
                                 n_maxit = 75)
  
  #format results
  results<-data.frame(permutation=i, 
                      model=c("univariate",
                              "multivariate",
                              "multivariate",
                              "multivariate",
                              "multivariate",
                              "multivariate"),
                      trait=c("SRCW",
                              "SRCW + Hard",
                              "SRCW + Dia",
                              "SRCW + QuadT",
                              "SRCW + QuadBC",
                              "SRCW	+ Hard + Dia + QuadT + QuadBC"),
                      r=c(one,
                          two,
                          three,
                          four,
                          five,
                          six))
  
  #bind in results
  cv_results<-rbind(cv_results, results)
  
  #remove
  remove(one,
         two,
         three,
         four,
         five,
         six,
         results)
  
  #announcement 
  print(paste("Permutation =", i, "complete"))
  print(Sys.time())
  
}

#save image
save.image(paste("UV-VS-MV-GS-", Sys.Date(), ".Rdata", sep = ""))

