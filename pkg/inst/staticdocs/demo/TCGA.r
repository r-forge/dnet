# This is a demo for TCGA mutational profile dataset from Kandoth et al
# 
# This dataset is available from TCGA (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/24132290" target="24132290">http://www.ncbi.nlm.nih.gov/pubmed/24132290</a>), containing somatic mutational profiles for 3096 cancer samples with survival data. These cancer samples belong to one of 12 major cancer types, including breast adenocarcinoma (BRCA), lung adenocarcinoma (LUAD), lung squamous cell carcinoma (LUSC), uterine corpus endometrial carcinoma (UCEC), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC), colon and rectal carcinoma (COAD/READ), bladder urothelial carcinoma (BLCA), kidney renal clear cell carcinoma (KIRC), ovarian serous carcinoma (OV) and acute myeloid leukaemia (LAML). For each patient sample, somatic mutations are represented as a profile of  states on genes, where non-zero entry indicates a gene for which how many mutations have occurred in the tumor relative to germ line. The dataset is provided as an 'ExpressionSet' object.
## assayData: exprs(TCGA_mutations), a matrix of 19171 genes X 3096 samples;
## phenoData: pData(TCGA_mutations), variables describing sample phenotypes (i.e. columns in assayData), including clinical/survival information about samples: "time" (i.e. survival time in days), "status" (i.e., survival status: 0=alive; 1=dead), "Age" (the patient age in years), "Gender" (the patient gender: male/female), "TCGA_tumor_type", "Tumor_stage", "Tumor_grade"
## featureData: fData(TCGA_mutations), variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol, "Desc" for gene description, "Synonyms" for gene symbol alias
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
list.pkg <- c("affy", "survival")
source("http://bioconductor.org/biocLite.R")
for(pkg in list.pkg){
    if(!require(pkg, character.only=T)){
        biocLite(pkg)
        lapply(pkg, library, character.only=T)
    }
}

# load an "ExpressionSet" object
load(url("http://dnet.r-forge.r-project.org/data/Datasets/TCGA_mutations.RData"))
#load("RData_Rd/data/Datasets/TCGA_mutations.RData")
eset <- TCGA_mutations
# extract information about phenotype data
pd <- pData(eset)
pd[1:5,]
# extract information about feature/gene data
fd <- fData(eset)
fd[1:5,]
# extract information about mutational data
md <- exprs(eset)
md[1:5,1:5]
# number of samples for each cancer type
table(pData(eset)$TCGA_tumor_type)
tumor_type <- sort(unique(pData(eset)$TCGA_tumor_type))

# Kaplan-Meier survival curves for individual tumor types
fit <- survfit(Surv(time, status) ~ TCGA_tumor_type, data=pd, type=c("kaplan-meier","fh")[1])
plot(fit, xscale=365.25, xlab = "Years", ylab="Survival", lty=1:2, col=rainbow(length(tumor_type))) 
legend("topright", tumor_type, lty=1:2, col=rainbow(length(tumor_type)))

# Survival analysis across tumor types using Cox proportional hazards model
# Cox regression yields an equation for the hazard/risk as a function of several explanatory variables
# Except for the gene mutational data in subject,  other explanatory variables (or called covariates) include: age, gender, and tumor type
## only those genes with mutations at least 1% of samples will be analysed
flag <- sapply(1:nrow(md), function(i) ifelse(sum(md[i,]!=0)>=0.01*ncol(md), T, F))
md_selected <- md[flag,]
## survival analysis to obtain hazard ratio (LR) and pvaules
gene_signif <- matrix(1, nrow=nrow(md_selected), ncol=2)
rownames(gene_signif) <- rownames(md_selected)
colnames(gene_signif) <- c("LR", "pvalue")
for(i in 1:nrow(md_selected)){
    if(i %% 100==0 || i==1 || i==nrow(md_selected)) message(sprintf("%d genes processed (%s)", i, as.character(Sys.time())), appendLF=T)
    ## fit a Cox proportional hazards model
    #gene_mut <- 1*(md_selected[i,]!=0)
    gene_mut <- md_selected[i,]
    data <- cbind(pd,gene=gene_mut)
    fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
    res <- as.matrix(anova(fit))
    ## 2nd: likelyhood ratio; 4th: pvalue
    gene_signif[i,] <- res[5,c(2,4)]
}
LR <- gene_signif[,1]
pvals <- gene_signif[,2]


# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 9.1). Only those associations with medium confidence (score>=400) are retained.
load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.string.RData"))
org.Hs.string
# restrict to those edges with high confidence (score>=700)
g <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=400])
g

# extract network that only contains genes in eset
ind <- match(V(g)$symbol, names(pvals))
## for extracted graph
nodes_mapped <- V(g)$name[!is.na(ind)]
network <- dNetInduce(g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
V(network)$name <- V(network)$symbol
#E(network)$weight <- E(network)$combined_score
network

# Identification of gene-active subnetwork
# 2) identification of gene-active subnetwork
#pvals <- pvals[V(network)$symbol]
## restrict the identified subnetwork to have the node size of 40 or so
#g <- dNetPipeline(g=network, pval=pvals, method="fdr", nsize=40)
g <- dNetPipeline(g=network, pval=pvals, method="fdr", fdr=3e-02)
g

B <- 100
gene_pvals_random <- matrix(1, nrow=nrow(md_selected), ncol=B)
rownames(gene_pvals_random) <- rownames(md_selected)
colnames(gene_pvals_random) <- 1:B
for(j in 1:B){
    for(i in 1:nrow(md_selected)){
        if(i %% 100==0 || i==1 || i==nrow(md_selected)) message(sprintf("Round %d: %d genes processed (%s)", j, i, as.character(Sys.time())), appendLF=T)
        # at least 1% of patients
        gene_mut <- md_selected[i,]
        ## fit a Cox proportional hazards model
        data <- cbind(pd,gene=gene_mut)
        ## do sampling (80% samples without replacement)
        flag <- sample(1:nrow(data), size=0.8*nrow(data), replace=F)
        data <- data[flag,]
        ## fit
        fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
        res <- as.matrix(anova(fit))
        ## 2nd: likelyhood ratio; 4th: pvalue
        gene_pvals_random[i,j] <- res[5,4]
    }
}

modules <- list()
for(j in 1:B){
    message(sprintf("###Jacknife resampling: %d replicate", j), appendLF=T)
    # For each jackknife replicate, recalculate survival p-values, which is used to calculate node scores. The same fdr as before is used
    #modules[[j]] <- dNetPipeline(g=network, pval=gene_pvals_random[,j], method="fdr", fdr=3e-02)
    modules[[j]] <- dNetPipeline(g=network, pval=gene_pvals_random[,j], method="fdr", nsize=vcount(g))
}

# append the confidence information from the source graphs into the target graph
cmodule <- dNetConfidence(target=g, sources=modules, plot=F)
# visualise the confidence target graph
visNet(cmodule, layout=glayout, edge.width=E(cmodule)$edgeConfidence/10, edge.label=E(cmodule)$edgeConfidence, edge.label.cex=0.6)





# 3) visualisation of the gene-active subnetwork itself
## the layout of the network visualisation (fixed in different visuals) 
glayout <- layout.fruchterman.reingold(g)
## color nodes according to communities (identified via a spin-glass model and simulated annealing)
com <- spinglass.community(g, spins=25)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(g, com)
## node sizes according to degrees
vdegrees <- igraph::degree(g)
## highlight different communities
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("grey", "black")[crossing(com,g)+1]
## visualise the subnetwrok
visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)


# fit a Cox proportional hazards model using the subnetwork
## for the whole network
data_g <- t(md[V(g)$name,])
data_g <- apply(data_g!=0, 1, sum)
data <- cbind(pd, net=data_g)
fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + net, data=data)
res <- as.matrix(anova(fit))
## 2nd: likelyhood ratio; 4th: pvalue
LR_g <- res[5,2]
pvals_g <- res[5,4]
## for the cumulative nodes from the network
cg_names <- names(sort(LR[V(g)$name], decreasing=T))
cg_signif <- matrix(1, nrow=length(cg_names), ncol=2)
rownames(cg_signif) <- cg_names
colnames(cg_signif) <- c("LR", "pvalue")
for(i in 1:length(cg_names)){
    data_g <- t(md[cg_names[1:i],])
    if(i!=1){
        data_g <- apply(data_g!=0, 1, sum)
    }else{
        data_g <- as.vector(data_g!=0)
    }
    data <- cbind(pd, cnet=data_g)
    fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + cnet, data=data)
    res <- as.matrix(anova(fit))
    cg_signif[i,] <- res[5,c(2,4)]
}
cg_signif[cg_signif[,2]==0,2] <- min(cg_signif[cg_signif[,2]!=0,2])


bp.LR.list <- list(All=LR, Neti=LR[cg_names], Netc=cg_signif[,1])
par(las=1, mar=c(5,8,4,2)) # all axis labels horizontal
boxplot(bp.LR.list, outline=F, horizontal=T, names=c("All genes", "Individual genes\n in the network", "Combined genes\n in the network"), col=c("red","green","blue"), border=par("fg"), xlab="Hazard Ratio")

plot(cg_signif[,1])
lines(LR[cg_names])

x <- 1:length(cg_names)
y <- LR[cg_names]
y <- cg_signif[,1]
df <- data.frame(x,y)
# model
mod <- lm(y ~ x, data = df)
# predicts + interval
newx <- df$x
preds <- predict(mod, newdata=data.frame(x=newx), interval='confidence')
# plot
plot(y ~ x, data=df, type='l')
# add fill
polygon(c(rev(newx),newx), c(rev(preds[,3]), preds[,2]), col='grey80', border=NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')



# 4) visualisation of the gene-active subnetwork overlaid by the node/gene score
max_colorbar <- ceiling(quantile(abs(V(g)$score),0.5))
visNet(g, glayout=glayout, pattern=V(g)$score, zlim=c(-1*max_colorbar,max_colorbar), vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 5) visualisation of the gene-active subnetwork overlaid by p-values
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- -1*log10(pvals[V(g)$name])
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="sphere", zlim=c(0,2), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 6) visualisation of the gene-active subnetwork overlaid by log-LR
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- log2(LR[V(g)$name])
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="sphere", zlim=c(-2,2), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 7) Network-based sample classifications and visualisations on 2D sample landscape
# it uses the gene-active subnetwork overlaid by all replication timing data
frac_mutated <- sapply(tumor_type, function(x) {
    e <- eset[, which(pData(eset)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)/ncol(e)
})
rownames(frac_mutated) <- fData(esetGene)$Symbol
frac_mutated[1:10,]

data <- frac_mutated[V(g)$name,]
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=3, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*3, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.8,border.color="888888", zlim=c(0,0.10), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 8) heatmap of replication timing data in the subnetwork
visHeatmapAdv(data, colormap=colormap, zlim=c(0,0.15), KeyValueName="Fraction of mutations")

# 9) output the subnetwork and their replication timing data
## Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(g)[,1], type="interaction", target=get.edgelist(g)[,2])
write.table(sif, file=paste("Survival_TCGA.sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
## Output the corresponding replication timing data
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste("Survival_TCGA.txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# 9) enrichment analysis for genes in the subnetwork
## get a list of genes in the subnetwork
data <- V(g)$name
data

## 9a) GOBP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="GOBP")
## visualise the top significant terms in the GOBP heirarchy
## first, load the GOBP ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOBP.RData"))
g <- ig.GOBP
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,]
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9b) GOMF enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="GOMF")
## visualise the top significant terms in the GOMF heirarchy
## first, load the GOMF ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
g <- ig.GOMF
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,]
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9c) MP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="MP", ontology.algorithm="elim")
## visualise the top significant terms in the MP heirarchy
## first, load the MP ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.MP.RData"))
g <- ig.MP
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$pvalue)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[nodes_query,2:3], cbind(zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])
## induce all possible paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info=c("none","term_id","term_name","both","full_term_name")[5], layout.orientation=c("left_right","top_bottom","bottom_top","right_left")[1], node.attrs=list(color=nodes.highlight))

## 9d) DO enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="DO", ontology.algorithm="pc")
## visualise the top significant terms in the DO heirarchy
## first, load the DO ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.DO.RData"))
g <- ig.DO
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[nodes_query,2:3], cbind(zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])
## induce all possible shortest paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,4), node.attrs=list(color=nodes.highlight))

## 9e) PS enrichment analysis
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="PS2")
## Loot at the evolution relevance along the path to the eukaryotic common ancestor
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), nOverlap=sapply(eTerm$overlap,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)

eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="SF")

eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2KEGG")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2CP")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2REACTOME")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC3TFT")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC3MIR")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC6")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC7")
nodes_query <- names(sort(eTerm$pvalue)[1:10])
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[nodes_query,2:3], cbind(nSet=sapply(eTerm$gs,length), nOverlap=sapply(eTerm$overlap,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])













# extract network that only contains genes in eset
ind <- match(V(g)$symbol, names(pvals))
## for extracted expression
esetGene <- eset[ind[!is.na(ind)],]
esetGene
## for extracted graph
nodes_mapped <- V(g)$name[!is.na(ind)]
network <- dNetInduce(g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=F)
V(network)$name <- V(network)$symbol
#E(network)$weight <- E(network)$combined_score
network

# For each gene, calculate number of samples (within a tumor type) having the mutated gene 
num_mutated <- sapply(tumor_type, function(x) {
    e <- esetGene[, which(pData(esetGene)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)
})
rownames(num_mutated) <- fData(esetGene)$Symbol
num_mutated[1:10,]

# For each gene, calculate fraction of samples (within a tumor type) having the mutated gene 
frac_mutated <- sapply(tumor_type, function(x) {
    e <- esetGene[, which(pData(esetGene)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)/ncol(e)
})
rownames(frac_mutated) <- fData(esetGene)$Symbol
frac_mutated[1:10,]

# define the "mutation sparseness" of genes in terms of a vector which stores the fraction of samples (within a tumor type) having the mutated gene 
# sparseness for a vector is: 1) one if the vector contains only a single non-zero value; 2) zero if and only if all elements are equal; 3) otherwise, the value interpolates smoothly between the two extremes
sparseness <- sapply(1:nrow(frac_mutated), function(i){
    v <- frac_mutated[i,]
    n <- length(v)
    norm1 <- sum(abs(v))
    norm2 <- sqrt(sum(v^2))
    (sqrt(n)-norm1/norm2) / (sqrt(n)-1)
})
sparseness <- matrix(sparseness, ncol=1)
rownames(sparseness) <- rownames(frac_mutated)
# derive the "mutation denseness" of genes
denseness <- 1- sparseness
hist(denseness,100)

# random walk with restart using denseness as seeds
PTmatrix <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=data.frame(frac_mutated, denseness), restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)

PTmatrix <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=frac_mutated, restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)


data <- t(PTmatrix)
n <- nrow(data)
res <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    res[i, ] <- (data[i, ] %*% t(data))
}
rownames(res) <- rownames(data)
colnames(res) <- rownames(data)
visHeatmapAdv(res, Rowv=F, Colv=F, zlim=c(0.0001,0.0004))

adjmatrix <- res > quantile(res[lower.tri(res)], 0.75)

sg <- graph.adjacency(adjmatrix, mode="undirected", weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
visNet(sg, edge.width=E(sg)$weight)



ind <- sample(1:nrow(frac_mutated))
seeds_random <- frac_mutated[ind,]
rownames(seeds_random) <- rownames(frac_mutated)
PT_random <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=seeds_random, restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)
data <- t(PT_random)
n <- nrow(data)
res_random <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    res_random[i, ] <- (data[i, ] %*% t(data))
}
rownames(res_random) <- rownames(data)
colnames(res_random) <- rownames(data)
res_random[lower.tri(res_random)]


a <- cbind(res[lower.tri(res)], res_random[lower.tri(res_random)])


# GSEA using PS
eTerm <- dGSEA(data=PTmatrix, identity="symbol", genome="Hs", ontology="PS", sizeRange=c(10,100000), which_distance=NULL, sigTail=c("two-tails","one-tail")[2])
which_sample=12
res <- dGSEAview(eTerm, which_sample=1, top_num=5, sortBy="nES", decreasing=T, details=TRUE)
visGSEA(eTerm, which_sample=which_sample, which_term=rownames(res)[1])
output <- dGSEAwrite(eTerm, which_content="nES", which_score="nES", cutoff=0, filename="eTerm.txt",keep.significance=F)

# visualise using advanced heatmap
data <- output[,6:ncol(output)]
rownames(data) <- paste(output$setID,output$setSize,output$name, sep="_")
visHeatmapAdv(data, Rowv=F, Colv=T, colormap="darkgreen-lightgreen-lightpink-darkred", zlim=c(0,2), margins = c(7,14),cexRow=1,cexCol=1)


