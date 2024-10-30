library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
WD<-getwd()
geneID2GO<-readMappings(file="geneID2GO")
geneUniverse <- names(geneID2GO)
genesOfInterest <- read.table("genesOfInterest",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1) 

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse


myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


myGOdata 

sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)

tutti<-usedGO(myGOdata)
raw_counts<-termStat(myGOdata, tutti)
write.table(raw_counts, file="raw_counts.txt")






resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")


mysummary <- summary(attributes(resultClassic)$score <= 0.01)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = numsignif)


output_file= "results"
sink(output_file)

output_file2 = paste(output_file,"Topgo", sep="_")
allRes
printGraph(myGOdata, resultParentchild, firstSigNodes = 20, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)


myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink() 



