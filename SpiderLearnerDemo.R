source("SpiderLearner.R")

library(affy)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedOvarianData")
library(curatedOvarianData)

standardize = function(x){return((x-mean(x))/sd(x))}

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)
lateStageClinicalInfo = pData(GSE32062.GPL6480_eset)

estrogen = read.table("~/Desktop/ovarianCarcinoma.txt",header=T,sep="\t")[,1]
lateStageSmall = lateStage[which(row.names(lateStage) %in% estrogen[-1]),] #GO_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY),]
lateStageSmall = t(lateStageSmall)
names(lateStageSmall) = colnames(lateStageSmall)
lateStageSmall = apply(lateStageSmall,2,standardize)

s = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new() 
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

candidates = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine)

for(candidate in candidates)
{
  s$addCandidate(candidate)
}

ptm <- proc.time()
slResults = s$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 1)
mod1time = proc.time() - ptm




