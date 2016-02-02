# Bioconductor routine for clustering together overlapping genes
library(GenomicFeatures)
txdb <- loadDb("genes.sqlite")
ebgene0 <- exonsBy(txdb, "gene")
select.genes <- unique(tx.db$GENEID[ idx ])
ebgene <- ebgene0[ select.genes ]
# => don't ignore strand for strand specific
fo <- findOverlaps(ebgene, ignore.strand=TRUE)
fo <- fo[queryHits(fo) < subjectHits(fo)]
mat <- as.matrix(fo)
library(graph)
library(RBGL)
myGraph <- ftM2graphNEL(mat, edgemode="undirected")
components <- connectedComp(myGraph)
components <- lapply(components, function(x) names(ebgene)[as.numeric(x)])
for (cluster in components) {
  tx.db$GENEID[tx.db$GENEID %in% cluster] <- paste0(cluster[1],"cluster")
}
