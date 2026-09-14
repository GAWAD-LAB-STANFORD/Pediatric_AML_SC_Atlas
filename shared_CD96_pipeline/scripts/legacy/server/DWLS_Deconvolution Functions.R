# ====================================================================
# DWLS_Deconvolution Functions.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript DWLS_Deconvolution Functions.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

#Deconvolution Functions
#install these packages if necessary...

if (!require("quadprog")) {
  install.packages("quadprog", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("reshape")) {
  install.packages("reshape", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("e1071")) {
  install.packages("e1071", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("Seurat")) {
  install.packages("Seurat", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("ROCR")) {
  install.packages("ROCR", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("varhandle")) {
  install.packages("varhandle", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("MAST")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("MAST")
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAST")
#load packages
library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)

#trim bulk and single-cell data to contain the same genes
trimData<-function(Signature,bulkData){
  Genes<-intersect(rownames(Signature),names(bulkData))
  B<-bulkData[Genes]
  S<-Signature[Genes,]
  return(list("sig"=S,"bulk"=B))
}


#solve using OLS, constrained such that cell type numbers>0
solveOLS<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#return cell number, not proportion
#do not print output
solveOLSInternal<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS<-function(S,B){
  #first solve OLS, use this solution to find a starting point for the weights
  solution<-solveOLSInternal(S,B)
  #now use dampened WLS, iterate weights until convergence
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-findDampeningConstant(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#solve WLS given a dampening constant
solveDampenedWLSj<-function(S,B,goldStandard,j){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[which(wsScaled>multiplier)]<-multiplier
  W<-diag(wsDampened)
  D<-t(S)%*%W%*%S
  d<- t(S)%*%W%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}

solveSVR<-function(S,B){
  #scaling
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1
  
  #perform SVR
  model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  print(round(coef/sum(coef),5))
  return(coef/sum(coef))
}

#perform DE analysis using Seurat
DEAnalysis<-function(Seurat_object,id,path){
  exprObj<-Seurat_object
  id_list<-as.list(unique(exprObj@meta.data$id))
  exprObj2<-SetIdent(exprObj,value=exprObj@meta.data$id)
  print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj2, ident.1 = i, ident.2 = NULL, 
                            only.pos = TRUE, test.use = "bimod")
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
  }
}

#build signature matrix using genes identified by DEAnalysis()
buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  
  #perform differential expression analysis
  DEAnalysis(scdata,id,path)
  
  numberofGenes<-c()
  for (i in unique(id)){
    load(file=paste(path,"/de_",i,".RData",sep=""))
    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),which(de_group$avg_logFC>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.",i,sep=""),de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }
  
  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
        Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1) #G is optimal gene number
  #
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
      Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  return(Sig)
}

##alternative differential expression method using MAST

#functions for DE

Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}

v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)
  
}  

#perform DE analysis using MAST	    
DEAnalysisMAST<-function(scdata,id,path){
  
  pseudo.count = 0.1
  data.used.log2   <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  for (i in unique(id)){
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))                          
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))   
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))
    
    DE <- bigtable[bigtable$log2_fc >diff.cutoff,] 
    dim(DE)
    if(dim(DE)[1]>1){
      data.1                 = data.used.log2[,cells.coord.list1]
      data.2                 = data.used.log2[,cells.coord.list2]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST 
      head(colData(vbeta.1))
      zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
      
      
      
      lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
      
      #. 4 save results
      write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
      save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))
    }
  }
}

#build signature matrix using genes identified by DEAnalysisMAST()
buildSignatureMatrixMAST<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  #compute differentially expressed genes for each cell type
  DEAnalysisMAST(scdata,id,path)
  
  #for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
  #is greater than 0.5
  numberofGenes<-c()
  for (i in unique(id)){
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
      load(file=paste(path,"/",i,"_MIST.RData", sep=""))
      pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr", n = length(cluster_lrTest.table[,3]))
      cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
      DEGenes<-cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted<pval.cutoff),which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
      nonMir = grep("MIR|Mir", DEGenes, invert = T)  # because Mir gene is usually not accurate 
      assign(paste("cluster_lrTest.table.",i,sep=""),cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
      numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
    }
  }
  
  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
        Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
      Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  return(Sig)
}

regroup.cor <- function(correlation.matrix, correlation.threshold = 0.95){
  
  top.cor <- which(abs(correlation.matrix) >= correlation.threshold & row(correlation.matrix) < col(correlation.matrix), arr.ind = TRUE)
  
  if(nrow(top.cor) != 0){
    
    ## reconstruct names from positions
    high_cor <- matrix(colnames(correlation.matrix)[top.cor], ncol = 2)
    
    count = 0; groups = list()
    
    for(elem in unique(high_cor[,1])){
      count = count + 1
      items = high_cor[high_cor[,1] %in% elem,]
      if(!is.null(dim(items))){items <- items %>% melt() %>% .[3] %>% .$value %>% unique()}
      groups[[paste("group",count,sep="_")]] <- items
    }
    
    for(elem in unique(high_cor[,2])){
      count = count + 1
      items = high_cor[high_cor[,2] %in% elem,]
      if(!is.null(dim(items))){items <- items %>% melt() %>% .[3] %>% .$value %>% unique()}
      groups[[paste("group",count,sep="_")]] <- items
    }
    
    # If one element present in multiple groups, keep largest
    groups <- groups[sapply(groups, length) %>% order() %>% rev()] #re-order NEEDED to ensure biggest groups and no repetition
    new.groups <- list()
    used <- c()
    count = 0
    for(i in names(groups)){
      count = count + 1
      to.regroup <- sapply(groups, function(x) groups[[i]] %in% x) %>% colSums() != 0
      to.regroup <- names(to.regroup)[to.regroup]
      if(sum(to.regroup %in% used) == 0){
        new.groups[[paste("group",count,sep="_")]] <- groups[to.regroup] %>% unlist() %>% unique()
        used <- c(used, to.regroup)
      }
      
    }
    
    annotation = melt(new.groups)
    rownames(annotation) <- annotation$value
    annotation$value <- NULL
    colnames(annotation) <- "new.group"
    
  } else {
    
    annotation = data.frame(new.group = "")
    
  }
  
  return(annotation)
  
}



collapsing <- function(input, duplicated.names){
  
  if(sum(duplicated(duplicated.names)) != 0){
    
    exp.var <- apply(input, 1, mean) %>% melt()
    exp.var$gene  <- duplicated.names
    exp.var$IQR <- apply(input, 1, IQR) %>% as.numeric()
    exp.var$row.number <- 1:nrow(exp.var)
    max.values = exp.var %>% 
      group_by(gene) %>% 
      summarise(across(c("value", "IQR"), ~ max(.x)))
    exp.var <- merge(max.values, exp.var, by = "gene")
    to.keep <- which(exp.var[,2] == exp.var[,4] & exp.var[,3] == exp.var[,5]) %>% sort()
    
    #if after this criteria, still duplicates, choose the first one (smallest row.number)!
    exp.var <- exp.var[to.keep,]
    min.row = exp.var %>% 
      group_by(gene) %>% 
      summarise(across(c("row.number"), ~ min(.x)))
    exp.var <- merge(min.row, exp.var, by = "gene")
    to.keep <- which(exp.var[,2] == exp.var[,7]) %>% sort()
    exp.var <- exp.var[to.keep,] %>% na.omit()
    input <- input[exp.var[,2], ]
    
    if(length(grep("ENS.*\\.[0-9]+", rownames(input))) > 1000){ # For datasets with ENSG names with format "ENSG00...003.14": keep genes with highest mean expression AND highest variability (first criterium insufficient, there are cases with same mean)
      
      rownames(input) <- strsplit(rownames(input),"\\.") %>% lapply(., function(x) x[1]) %>% unlist()
      
    } else if(length(grep("ENS.*\\__", rownames(input))) > 1000){  # For datasets with gene names as "ENSG00...__HGNC": remove the HGNC part
      
      rownames(input) <- strsplit(rownames(input),"__") %>% lapply(., function(x) x[1]) %>% unlist()
      
    } else {
      
      rownames(input) <- exp.var$gene
      
    }
    
  }
  
  return(input)
  
}



transformation2 <- function(X, Y, leave.one.out = TRUE) {
  
  # use the same genes for all input datasets
  Genes <- intersect(row.names(Y), row.names(X))
  
  X <- as.matrix(X[Genes,])
  Y <- as.matrix(Y[Genes,])
  
  if(leave.one.out){
    
    pred <- intersect(colnames(X), colnames(Y)) # matching samples with sc and bulk data
    X.new <- matrix(0, nrow=dim(X)[1], ncol=length(pred))
    
    for(j in 1:length(pred)){
      
      X.train <- as.matrix(X[,pred[-j]])
      X.test <- as.matrix(X[,pred[j]])
      Y.train <- as.matrix(Y[,pred])
      
      # track the mean and sd after leaving one out:
      X.mean <- rowMeans(X.train)
      X.sd <- apply(X.train,1,sd) 
      Y.mean <- rowMeans(Y)
      Y.sd <- apply(Y,1,sd)
      
      for (i in 1:dim(X)[1]){
        
        # transforming over genes
        x <- X.test[i]
        sigma_j <- Y.sd[i]*sqrt((length(Y[i,])-1)/(length(Y[i,])+1))
        x.new <- (x-X.mean[i])/X.sd[i]
        x.new <- x.new*sigma_j + Y.mean[i]
        
        X.new[i,j] <- x.new
        
      }
      
    }
    
  } else {
    
    # transforming over genes
    X.new <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
    
    for (i in 1:dim(X)[1]){
      
      y <- Y[i,]
      x <- X[i,]
      l <- length(y)
      sigma_j <- sd(y)*sqrt((l-1)/(l+1))
      x.new <- (x-mean(x))/sd(x)
      x.new <- x.new*sigma_j + mean(y)
      X.new[i,] <- x.new
      
    }
    
  }
  
  # explicit non-negativity constraint:
  X.new = apply(X.new,2,function(x) ifelse(x < 0, 0, x))
  
  rownames(X.new) = rownames(X); colnames(X.new) = colnames(X)
  
  return(X.new)
  
}



bulkC.fromSC <- function(scC, phenoDataC){
  
  phenoDataC = phenoDataC[colnames(scC),]
  cellType <- phenoDataC$cellType
  group = list()
  
  for(i in unique(cellType)){ 
    
    group[[i]] <- which(cellType %in% i)
    
  }
  
  C = lapply(group,function(x) Matrix::rowMeans(scC[,x])) 
  C = do.call(cbind.data.frame, C)
  
  return(as.matrix(C))
  
}



Scaling <- function(matrix, option, phenoDataC=NULL){
  
  #Avoid Error: Input matrix x contains at least one null or NA-filled row.
  matrix = matrix[rowSums(matrix) != 0,]
  #Avoid error if all elements within a row are equal (e.g. all 0, or all a common value after log/sqrt/vst transformation)
  matrix = matrix[!apply(matrix, 1, function(x) var(x) == 0),]
  matrix = matrix[,colSums(matrix) != 0]
  
  if (option == "LogNormalize"){
    
    matrix = expm1(Seurat::LogNormalize(matrix, verbose = FALSE)) %>% as.matrix()
    
  } else if (option == "TMM"){# CPM counts coming from TMM-normalized library sizes; https://support.bioconductor.org/p/114798/
    
    if(!is.null(phenoDataC)){
      
      Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])
      
    } else {
      
      Celltype = colnames(matrix)
      
    }
    
    matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
    CtrlGenes <- grep("ERCC-",rownames(data))
    
    if(length(CtrlGenes) > 1){
      
      spikes <- data[CtrlGenes,]
      spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
      matrix$samples$norm.factors <- spikes$samples$norm.factors
      
    } else {
      
      matrix <- edgeR::calcNormFactors(matrix, method = "TMM")  
      
    }
    
    matrix <- edgeR::cpm(matrix)
    
  } else if (option == "TPM"){
    
    # MULTI-CORE VERSION: adapted from #devtools::source_url('https://raw.githubusercontent.com/dviraran/SingleR/master/R/HelperFunctions.R', sha1 = "df5560b4ebb28295349aadcbe5b0dc77d847fd9e")
    
    TPM <- function(counts,lengths=NULL){
      
      require(foreach)
      require(Matrix)
      
      if (is.null(lengths)) {
        data('gene_lengths')
      }
      
      A = intersect(rownames(counts),names(lengths))
      counts = as.matrix(counts[A,])
      lengths = lengths[A]
      rate = counts / lengths 
      
      num_cores <- min(4, parallel::detectCores()) #slurm will automatically restrict this to parameter --cpus-per-task
      print(paste("num_cores = ", num_cores, sep = ""))
      
      doMC::registerDoMC(cores = num_cores)
      jump = 100
      tpm <- foreach::foreach(elem = seq(1,ncol(counts), jump), .combine='cbind.data.frame') %dopar% {   
        
        lsizes <- colSums(rate[, elem:min(ncol(counts),(elem + jump - 1))])
        lo <- 1e6*rate[, elem:min(ncol(counts),(elem + jump - 1))]
        tpm <- lo %*% diag(1/lsizes)
        tpm
        
      }
      
      colnames(tpm) <- colnames(counts)
      return(tpm)
      
    }
    
    if(! file.exists("human_lengths.rda")){
      download.file("https://github.com/dviraran/SingleR/blob/master/data/human_lengths.rda?raw=true","./human_lengths.rda")
    }
    load("./human_lengths.rda")
    
    # Doesn't work with Ensembl IDs:
    if(length(grep("ENSG000",rownames(matrix))) > 50){
      
      suppressMessages(library("AnnotationDbi"))
      suppressMessages(library("org.Hs.eg.db"))
      temp = mapIds(org.Hs.eg.db,keys = names(human_lengths), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
      names(human_lengths) = as.character(temp)
      
    }
    
    matrix = TPM(counts = matrix, lengths = human_lengths)
    rownames(matrix) = toupper(rownames(matrix))
    
  } else if (option == "TPM.murine"){
    
    # MULTI-CORE VERSION: adapted from #devtools::source_url('https://raw.githubusercontent.com/dviraran/SingleR/master/R/HelperFunctions.R', sha1 = "df5560b4ebb28295349aadcbe5b0dc77d847fd9e")
    
    TPM <- function(counts,lengths=NULL){
      
      require(foreach)
      require(Matrix)
      
      if (is.null(lengths)) {
        data('gene_lengths')
      }
      
      
      A = intersect(rownames(counts),names(lengths))
      counts = as.matrix(counts[A,])
      lengths = lengths[A]
      rate = counts / lengths 
      
      num_cores <- min(4, parallel::detectCores()) #slurm will automatically restrict this to parameter --cpus-per-task
      print(paste("num_cores = ", num_cores, sep = ""))
      
      doMC::registerDoMC(cores = num_cores)
      jump = 100
      tpm <- foreach::foreach(elem = seq(1,ncol(counts), jump), .combine='cbind.data.frame') %dopar% {   
        
        lsizes <- colSums(rate[, elem:min(ncol(counts),(elem + jump - 1))])
        lo <- 1e6*rate[, elem:min(ncol(counts),(elem + jump - 1))]
        tpm <- lo %*% diag(1/lsizes)
        tpm
        
      }
      
      colnames(tpm) <- colnames(counts)
      return(tpm)
      
    }
    
    if(! file.exists("mouse_lengths.rda")){
      download.file("https://github.com/dviraran/SingleR/blob/master/data/mouse_lengths.rda?raw=true","./mouse_lengths.rda")
    }
    
    load("./mouse_lengths.rda")
    
    ## Append version with ENSMUSG000 names
    mouse_lengths2 <- mouse_lengths
    
    # Doesn't work with Ensembl IDs:
    suppressMessages(library("AnnotationDbi"))
    suppressMessages(library("org.Mm.eg.db"))
    temp = mapIds(org.Mm.eg.db,keys = names(mouse_lengths), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
    names(mouse_lengths) = as.character(temp)
    
    mouse_lengths = c(mouse_lengths,mouse_lengths2)
    matrix = TPM(counts = matrix, lengths = mouse_lengths)
    
    ####################################################################################
    ## scRNA-seq specific  
    
  } else if (option == "SCTransform"){
    
    matrix = as(matrix, "dgCMatrix")
    matrix = sctransform::vst(matrix, return_corrected_umi = TRUE, show_progress = FALSE)$umi_corrected
    
  } else if (option == "scran"){
    
    sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(matrix)))
    sce = scran::computeSumFactors(sce, clusters=NULL)
    sce = scater::logNormCounts(sce, log = FALSE)
    matrix = normcounts(sce)
    
  } else if (option == "scater"){  
    
    sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(matrix)))
    size_factors = scater::librarySizeFactors(sce)
    sce = scater::logNormCounts(sce, log = FALSE)
    matrix = normcounts(sce)
    
  }
  
  return(matrix)
  
}



Deconvolution <- function(T, C, method, phenoDataC, P = NULL, elem = NULL, STRING = NULL, marker_distrib, refProfiles.var){ 
  
  bulk_methods = c("CIBERSORT","nnls","FARDEEP","RLR")
  sc_methods = c("MuSiC","MuSiC_with_markers","DWLS","SQUID", "Bisque", "Bisque_with_markers")
  
  ########## Using marker information for bulk_methods
  
  if(method %in% bulk_methods){
    
    C = C[rownames(C) %in% unique(marker_distrib$gene), , drop = FALSE]
    T = T[rownames(T) %in% unique(marker_distrib$gene), , drop = FALSE]
    refProfiles.var = refProfiles.var[rownames(refProfiles.var) %in% unique(marker_distrib$gene), , drop = FALSE]
    
  } else { ### For scRNA-seq methods 
    
    if(length(grep("[N-n]ame",colnames(phenoDataC))) > 0){
      sample_column = grep("[N-n]ame",colnames(phenoDataC))
    } else {
      sample_column = grep("[S-s]ample|[S-s]ubject",colnames(phenoDataC))
    }
    
    colnames(phenoDataC)[sample_column] = "SubjectName"
    rownames(phenoDataC) = phenoDataC$cellID
    # establish same order in (sc)C and phenoDataC:
    phenoDataC <- phenoDataC[match(colnames(C),phenoDataC$cellID),]
    
    if(method %in% c("MuSiC","MuSiC_with_markers", "Bisque", "Bisque_with_markers")){
      
      require(xbioc)
      C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C), phenoData = Biobase::AnnotatedDataFrame(phenoDataC))
      T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T))
      
    }
    
  }
  
  ##########    MATRIX DIMENSION APPROPRIATENESS    ##########
  
  keep = intersect(rownames(C),rownames(T)) 
  C = C[keep, , drop = FALSE]
  T = T[keep, , drop = FALSE]
  
  ###################################
  
  if(method == "CIBERSORT"){ 
    
    #source("./CIBERSORT.R")
    RESULTS = CIBERSORT(sig_matrix = C, mixture_file = T, QN = FALSE) 
    RESULTS = t(RESULTS[,1:(ncol(RESULTS)-3),drop=FALSE]) 
    
  } else if (method == "nnls"){
    
    require(nnls)
    RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- colnames(C)
    
  } else if (method == "FARDEEP"){
    
    require(FARDEEP)
    RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    # These 2 lines are similar as retrieving $relative.beta instead of abs.beta + re-scaling
    
  } else if (method == "RLR"){ #RLR = robust linear regression
    
    require(MASS)
    RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
    
    ###################################
    ###################################
    
  } else if (method == "MuSiC"){
    
    require(MuSiC)
    RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                  markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                                  verbose = F)$Est.prop.weighted)
    
  } else if (method == "MuSiC_with_markers"){
    
    require(MuSiC)
    RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                  markers = unique(marker_distrib$gene), normalize = FALSE, samples = 'SubjectName', 
                                  verbose = F)$Est.prop.weighted)
    
  } else if (method == "Bisque"){#By default, Bisque uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter
    
    require(BisqueRNA)
    RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers = NULL, use.overlap = FALSE)$bulk.props 
    
  } else if (method == "Bisque_with_markers"){#By default, Bisque uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter
    
    require(BisqueRNA)
    RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers = marker_distrib$gene, use.overlap = FALSE)$bulk.props 
    
  } else if (method == "DWLS"){
    
    require(DWLS)
    path = paste(getwd(),"/results_",STRING,sep="")
    
    if(! dir.exists(path)){ dir.create(path) } #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
    
    if(!file.exists(paste(path,"Sig.RData",sep="/"))){
      
      Signature <- DWLS::buildSignatureMatrixMAST(scdata = C, id = as.character(phenoDataC$cellType), path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)
      
    } else {#re-load signature and remove CT column + its correspondent markers
      
      load(paste(path,"Sig.RData",sep="/"))
      Signature <- Sig
      
    }
    
    Signature = as.matrix(Signature)
    
    RESULTS <- apply(T, 2, function(x){
      b = setNames(x, rownames(T))
      tr <- DWLS::trimData(Signature, b)
      RES <- t(DWLS::solveDampenedWLS(tr$sig, tr$bulk))
    })
    
    rownames(RESULTS) <- as.character(unique(phenoDataC$cellType))
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    print(head(RESULTS))
    
  } else if (method == "SQUID"){ 
    
    # Transforming T :
    Z = bulkC.fromSC(scC = C, phenoDataC = phenoDataC)
    X = T
    P <- as.matrix(P[base::colnames(Z),,drop = FALSE])
    Y = Z %*% P
    
    X.new = transformation2(X = X, Y = Y, leave.one.out = FALSE)
    X.new[!is.finite(X.new)] <- 0 
    
    # take common genes
    Genes <- intersect(rownames(Z),rownames(X.new))
    
    X.new <- as.matrix.Vector(X.new[Genes,])
    Z <- as.matrix.Vector(Z[Genes,])
    
    RESULTS <- apply(X.new, 2, function(x){
      RES <- t(DWLS::solveDampenedWLS(S = Z, B = x))
    })
    
    rownames(RESULTS) <- colnames(Z)
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  }
  
  RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)), , drop = FALSE] 
  RESULTS = data.table::melt(RESULTS)
  colnames(RESULTS) <- c("cell_type","tissue","observed_fraction")
  RESULTS$cell_type = as.character(RESULTS$cell_type)
  RESULTS$tissue = as.character(RESULTS$tissue)
  
  if(!is.null(P)){
    
    P = P[gtools::mixedsort(rownames(P)),,drop = FALSE] %>% data.frame(., check.names = FALSE)
    P$cell_type = rownames(P)
    P = data.table::melt(P, id.vars="cell_type")
    colnames(P) <-c("cell_type","tissue","expected_fraction")
    P$cell_type = as.character(P$cell_type)
    P$tissue = as.character(P$tissue)
    
    RESULTS = merge(RESULTS, P, by = c("cell_type", "tissue"), all = TRUE)
    RESULTS[is.na(RESULTS)] <- 0
    RESULTS$expected_fraction <- round(RESULTS$expected_fraction, 3)
    RESULTS$observed_fraction <- round(RESULTS$observed_fraction, 3)
    
  }
  
  return(RESULTS) 
  
}

