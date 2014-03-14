# This is a demo for human embryo dataset from Fang et al
# 
# This human embryo expression dataset (available from <a href="http://www.ncbi.nlm.nih.gov/pubmed/20643359" target="20643359">http://www.ncbi.nlm.nih.gov/pubmed/20643359</a>) involves six successive developmental stages (S9-S14) with three replicates (R1-R3) for each stage, including:
## Fang: an expression matrix of 5,441 genes X 18 samples;
## Fang.geneinfo: a matrix of 5,441 X 3 containing gene information;
## Fang.sampleinfo: a matrix of 18 X 3 containing sample information.
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
if(!require(affy)){
    install.packages("affy",repos="http://www.stats.bris.ac.uk/R")
    library(affy)
}
if(!require(limma)){
    install.packages("limma",repos="http://www.stats.bris.ac.uk/R")
    library(limma)
}

# import data containing three variables ('Fang', 'Fang.geneinfo' and 'Fang.sampleinfo')
data(Fang)
data <- Fang
fdata <- as.data.frame(Fang.geneinfo[,2:3])
rownames(fdata) <- Fang.geneinfo[,1]
pdata <- as.data.frame(Fang.sampleinfo[,2:3])
rownames(pdata) <- Fang.sampleinfo[,1]

# create the 'eset' object
colmatch <- match(rownames(pdata),colnames(data))
rowmatch <- match(rownames(fdata),rownames(data))
data <- data[rowmatch,colmatch]
eset <- new("ExpressionSet", exprs=as.matrix(data), phenoData=as(pdata,"AnnotatedDataFrame"), featureData=as(fdata,"AnnotatedDataFrame"))

# A function to convert probeset-centric to entrezgene-centric expression levels
prob2gene <- function(eset){
    fdat <- fData(eset)
    tmp <- as.matrix(unique(fdat[c("EntrezGene", "Symbol")]))
    forder <- tmp[order(as.numeric(tmp[,1])),]
    forder <- forder[!is.na(forder[,1]),]
    rownames(forder) <- forder[,2]
    system.time({
        dat <- exprs(eset)
        edat <- matrix(data=NA, nrow=nrow(forder), ncol=ncol(dat))
        for (i in 1:nrow(forder)){
            ind <- which(fdat$EntrezGene==as.numeric(forder[i,1]))
            if (length(ind) == 1){
                edat[i,] <- dat[ind,]
            }else{
                edat[i,] <- apply(dat[ind,],2,mean)
            }
        }
    })
    
    rownames(edat) <- rownames(forder) # as gene symbols
    colnames(edat) <- rownames(pData(eset))
    esetGene <- new('ExpressionSet',exprs=data.frame(edat),phenoData=as(pData(eset),"AnnotatedDataFrame"),featureData=as(data.frame(forder),"AnnotatedDataFrame"))
    return(esetGene)
}
esetGene <- prob2gene(eset)
esetGene

# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 9.0.5). Only those associations with medium confidence (score>=0.4) are retained.
load(url("http://dnet.r-forge.r-project.org/data/org.Hs.string.RData"))
org.Hs.string

# extract network that only contains genes in esetGene
ind <- match(V(org.Hs.string)$symbol, rownames(esetGene))
## for extracted expression
esetGeneSub <- esetGene[ind[!is.na(ind)],]
esetGeneSub
## for extracted graph
nodes_mapped <- V(org.Hs.string)$name[!is.na(ind)]
network <- dNetInduce(g=org.Hs.string, nodes_query=nodes_mapped, knn=0, remove.loops=T, largest.comp=T)
V(network)$name <- V(network)$symbol
network

# 1) preparation of node p-values
## define the design matrix in a order manner
all <- as.vector(pData(esetGene)$Stage)
level <- levels(factor(all))
index_level <- sapply(level, function(x) which(all==x)[1])
level_sorted <- all[sort(index_level, decreasing=F)]
design <- sapply(level_sorted, function(x) as.numeric(all==x)) # Convert a factor column to multiple boolean columns
## a linear model is fitted for every gene by the function lmFit
fit <- lmFit(exprs(esetGene), design)
## define a contrast matrix
contrasts <- dContrast(level_sorted, contrast.type="average")
contrast.matrix <- makeContrasts(contrasts=contrasts$each, levels=design)
colnames(contrast.matrix) <- contrasts$name
## computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# for p-value
pvals <- as.matrix(fit2$p.value)
# for adjusted p-value
adjpvals <- sapply(1:ncol(pvals),function(x) {
    p.adjust(pvals[,x], method="BH")
})
colnames(adjpvals) <- colnames(pvals)
# num of differentially expressed genes
apply(adjpvals<1e-2, 2, sum)

# only for the comparisons of S9 against the average
my_contrast <- "S9"
# get the p-values and calculate the scores thereupon
pval <- pvals[,my_contrast]
#pval <- dPvalAggregate(pvals, method="orderStatistic", order=ncol(pvals))
#pval <- dPvalAggregate(pvals, method="fishers", order=ncol(pvals))

data <- exprs(esetGeneSub)[adjpvals[,my_contrast]<1e-2,]
heatmap(as.matrix(data),col=visColormap("bwr")(64),zlim=c(min(data),max(data)), scale="none", cexRow=0.2+0.5/log10(nrow(data)), cexCol=0.2+0.5/log10(ncol(data)), Rowv=NULL,Colv=NA)

# 2) identification of module
g <- dNetPipeline(g=network, pval=pval, nsize=40)
glayout <- layout.fruchterman.reingold(g)

# 3) color nodes according to communities identified via a spin-glass model and simulated annealing
#com <- spinglass.community(g, spins=3)
com <- walktrap.community(g, modularity=T)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(g, com)

# 4) size nodes according to degrees
vdegrees <- igraph::degree(g)

# 5) sort nodes: first by communities and then degrees
tmp <- data.frame(ind=1:vcount(g), vgroups, vdegrees)
ordering <- tmp[order(vgroups,vdegrees),]$ind

# 6) visualise graph using 1-dimensional arc diagram
visNetArc(g, ordering=ordering, labels=V(g)$geneSymbol, vertex.label.color=vcolors, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.size=log(vdegrees)+0.1, vertex.label.cex=0.4)

# 7) visualise graph using circle diagram
# 7a) drawn into a single circle 
visNetCircle(g=g, com=com, ordering=ordering, colormap=colormap, vertex.label=V(g)$symbol, vertex.size=igraph::degree(g)+5, vertex.label.color="black", vertex.label.cex=0.6, vertex.label.dist=0.75, vertex.shape="sphere", edge.color.within="grey", edge.color.crossing="black", edge.width=1, edge.lty=1, mark.shape=1, mark.expand=10)
# 7b) drawn into multiple circles
visNetCircle(g=g, com=com, circles="multiple", ordering=ordering, colormap=colormap, vertex.label=V(g)$symbol, vertex.size=igraph::degree(g)+5, vertex.label.color="black", vertex.label.cex=0.6, vertex.label.dist=0.25, vertex.shape="sphere", edge.color.within="grey", edge.color.crossing="black", edge.width=1, edge.lty=1, mark.shape=1, mark.expand=10)

# 8) as comparison, also visualise graph on 2-dimensional layout 
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("grey", "black")[crossing(com,g)+1]
visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

legend_name <- paste("C",1:length(mcolors)," (n=",com$csize,", pval=",signif(com$significance,digits=2),")",sep='')
legend("bottomleft", legend=legend_name, fill=mcolors, bty="n", cex=0.6)

# 9) color by score and FC
# colored by score
visNet(g, glayout=glayout, pattern=V(g)$score, zlim=c(-1*ceiling(max(abs(V(g)$score))),ceiling(max(abs(V(g)$score)))), vertex.shape="circle", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)
# colored by FC
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- fit2$coefficients[V(g)$name,my_contrast]
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 10) color by additional data
##data <- exprs(esetGene)[V(g)$name,]
data <- as.matrix(fit2$coefficients[V(g)$name,])
visNetMul(g=g, data=data, height=ceiling(sqrt(ncol(data)))*2, newpage=T,glayout=glayout,colormap=colormap,vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.8,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 11) color by additional data (be reordered)
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=2, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*2, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.8,border.color="888888", mark.groups=mark.groups, mark.col=mark.col, mark.border=NA, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 12) heatmap of subnetwork
data[data < -2] <- -2
data[data > 2] <- 2
heatmap(as.matrix(data),col=visColormap(colormap=colormap)(64),zlim=c(-2,2), scale="none", cexRow=0.2+0.5/log10(nrow(data)), cexCol=0.2+0.5/log10(ncol(data)), Rowv=NULL,Colv=NULL)
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste(my_contrast,".txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# 13) Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(g)[,1], type="interaction", target=get.edgelist(g)[,2])
write.table(sif, file=paste(my_contrast,".sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
