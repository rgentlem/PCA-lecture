
library(scater)
load("sc_exp_norm_trans.rda")
rotation = attr(reducedDims(sc_exp_norm_trans)$PCA, "rotation")
gn = row.names(rotation)

library(org.Hs.eg.db)
library("Matrix")
##get all the GO mappings for the gene names used in our PCA - gn
ss1 = select(org.Hs.eg.db, keys=gn, columns=c("SYMBOL","GOALL"), keytype="SYMBOL")
ss2 = split(ss1$GOALL, ss1$SYMBOL)

##map from GO to symbol - so we can regress
## we will require 10 genes for any GO term to be included
## probably we should eliminate nodes with too many genes as well...
##
##ss3 is a list, each element of the list is the names/symbols
## of the gene that are attached to the GO category
ss3 = split(ss1$SYMBOL, ss1$GOALL)
ss3len = sapply(ss3, length)


ss4 = ss3[ss3len>9 & ss3len<50]
##FIXME - how to create a sparse matrix from ss4
## we are going to have columns as GO terms, rows as Genes
##FIXME - why are there multiple matches per GO term - below
numcols = length(ss4)
colNames = names(ss4)
rowLabs = unique(unlist(ss4))
p= vector("list", length=numcols)
names(p) = names(ss4)
for(i in 1:numcols)
  p[[i]] = match(unique(ss4[[i]]), rowLabs)

GenesByGO = sparseMatrix(j=rep(1:numcols, sapply(p, length)), i=unlist(p),x=1,
                         dimnames=list(Genes=rowLabs,GO=colNames))
GBG = as.matrix(GenesByGO)
save(GBG, file="GBG.rda")
save(ss4, file="ss4.rda")