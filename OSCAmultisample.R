##based on this vignette
##https://bioconductor.org/packages/devel/bioc/vignettes/batchelor/inst/doc/correction.html


library(scRNAseq)
sce1 <- ZeiselBrainData()
sce1

sce2 <- TasicBrainData()
sce2

library(scuttle)
library(batchelor)

sce1 <- addPerCellQC(sce1, subsets=list(Mito=grep("mt-", rownames(sce1))))
qc1 <- quickPerCellQC(colData(sce1), sub.fields="subsets_Mito_percent")
sce1 <- sce1[,!qc1$discard]

sce2 <- addPerCellQC(sce2, subsets=list(Mito=grep("mt_", rownames(sce2))))
qc2 <- quickPerCellQC(colData(sce2), sub.fields="subsets_Mito_percent")
sce2 <- sce2[,!qc2$discard]

universe <- intersect(rownames(sce1), rownames(sce2))
sce1 <- sce1[universe,]
sce2 <- sce2[universe,]

out <- multiBatchNorm(sce1, sce2)
sce1 <- out[[1]]
sce2 <- out[[2]]

library(scran)
dec1 <- modelGeneVar(sce1)
dec2 <- modelGeneVar(sce2)
combined.dec <- combineVar(dec1, dec2)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

combined <- correctExperiments(A=sce1, B=sce2, PARAM=NoCorrectParam())

library(scater)
#set.seed(100)
#combined <- runPCA(combined, subset_row=chosen.hvgs)
#combined <- runTSNE(combined, dimred="PCA")
#plotTSNE(combined, colour_by="batch")

mergedE = assays(combined)[[1]]
mergedE = mergedE[chosen.hvgs,]

pcaComb = prcomp(t(mergedE), scale=TRUE)

##look at first 9 PCs
multR2raw = sapply(1:9, function(x) {
  summary(lm(pcaComb$x[,x]~ combined$batch - 1))$adj.r.squared
})
round(multR2raw, digits=4)

##from the regression we see that the batch effects are pretty substantial
##and from their tsne plots

#pca_results <- reducedDims(combined)$PCA
#for(i in 1:20) print(summary(lm(pca_results[,i] ~ combined$batch))$r.squared)

## now for their batch correction

set.seed(101)
f.out <- fastMNN(A=sce1, B=sce2, subset.row=chosen.hvgs)

pcc = prcomp(t(assays(f.out)[[1]]), scale=T)

##and things look better - but this is a complex experiment - with different types of cells and 
## that means a deeper look is going to be needed
multR2mnn = sapply(1:9, function(x) {
  summary(lm(pcc$x[,x]~ f.out$batch - 1))$adj.r.squared
})
round(multR2mnn, digits=4)

##this is probably the output of multiBatchPCA
#pcar2 = reducedDims(f.out)$corrected

save(pcaComb, file="pcaComb.rda")
save(pcc, file="pcc.rda")
batch=f.out$batch
save(batch, file="batch.rda")



