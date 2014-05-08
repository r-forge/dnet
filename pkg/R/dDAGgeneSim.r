#' Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dDAGgeneSim} is supposed to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param genes the genes between which pair-wise semantic similarity is calculated. If NULL, all genes annotatable in the input dag will be used for calcluation, which is very prohibitively expensive!
#' @param method.gene the method used for how to derive semantic similarity between genes from semantic similarity between terms. It can be "average" for average similarity between any two terms (one from gene 1, the other from gene 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each gene in a gene pair, first calculate maximum similarity for each term, then their average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two genes in pair), "BM.max" for BM based maximum similarity (i.e. for each gene in a gene pair, first calculate maximum similarity for each term, then their average of maximum similarity; the final BM-based maximum similiary is the maximum of the pre-calculated average between two genes in pair), "BM.complete" for BM based complete-linkage similarity (inspired by complete-linkage procedure: for each gene in a gene pair, first calculate maximum similarity for each term; then take the mimumum of pre-calculated maximum similarity between two genes in pair). By default, it sets "BM.average"
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative information ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{http://webdocs.cs.ualberta.ca/~lindek/papers/sim.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - diference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param force logical to indicate whether the only most specific terms (for each gene) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{dDAGannotate}}).
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix containing pair-wise semantic similarity between input terms
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @import Matrix
#' @seealso \code{\link{dDAGtermSim}}, \code{\link{dDAGinduce}}, \code{\link{dDAGtip}}
#' @include dDAGgeneSim.r
#' @examples
#' \dontrun{
#' # 1) load HPPA as igraph object
#' data(ig.HPPA)
#' g <- ig.HPPA
#'
#' # 2) load human genes annotated by HPPA
#' data(org.Hs.egHPPA)
#'
#' # 3) prepare for ontology and its annotation information 
#' dag <- dDAGannotate(g, annotations=org.Hs.egHPPA, path.mode="all_paths", verbose=TRUE)
#'
#' # 4) calculate pair-wise semantic similarity between 5 randomly chosen terms 
#' allgenes <- unique(unlist(V(dag)$annotations))
#' genes <- sample(allgenes,5)
#' sim <- dDAGgeneSim(g=dag, genes=genes, method.gene=c("average","max","BM.average","BM.max","BM.complete")[3], method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita")[3], verbose=TRUE)
#' sim
#' }

dDAGgeneSim <- function (g, genes=NULL, method.gene=c("average","max","BM.average","BM.max","BM.complete")[3], method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita")[3], force=TRUE, verbose=TRUE)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################

    method.gene <- match.arg(method.gene)
    method.term <- match.arg(method.term)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(is.null(V(ig)$annotations) | is.null(V(ig)$IC)){
        stop("The function requires that input graph has already contained annotation data. Please first run 'dDAGannotate'.\n")
    }

    ####################################################
    # A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    ####################################################

    if(verbose){
        message(sprintf("First, extract all annotatable genes (%s)...", as.character(Sys.time())), appendLF=T)
    }

    anno <- V(ig)$annotations
    allgenes <- sort(as.numeric(unique(unlist(anno))))
    
    ## checking input genes
    genes <- genes[!is.na(genes)]
    if(is.null(genes) || is.na(genes)){
        genes <- allgenes
    }else{
        flag <- genes %in% allgenes
        if(sum(flag)!=0){
            genes <- genes[flag]
        }else{
            genes <- allgenes
        }
    }
    
    ## pre-compute a sparse matrix of input genes x terms
    allterms <- 1:length(anno)
    sGT <- Matrix::Matrix(0, nrow=length(genes), ncol=length(allterms), sparse=T)
    for(j in 1:length(allterms)){
        ind <- match(anno[[j]], genes)
        flag <- ind[!is.na(ind)]
        if(length(flag)!=0){
            sGT[flag,j] <- 1
        }
    }
    colnames(sGT) <- V(ig)$name
    rownames(sGT) <- genes

    if(verbose){
        message(sprintf("\tthere are %d input genes amongst %d annotatable genes", length(genes), length(allgenes)), appendLF=T)
    }
    
    ## a list of genes, each containing terms annotated by
    genes2terms <- sapply(1:length(genes), function(x){
        res <- names(which(sGT[x,]==1))
        if(force){
            subg <- dDAGinduce(ig, nodes_query=res, path.mode="all_paths")
            res <- dDAGtip(subg)
        }
        return(res)
    })
    names(genes2terms) <- genes
    terms <- unique(unlist(genes2terms))
    
    if(verbose){
        if(force){
            message(sprintf("Second, pre-compute semantic similarity between %d terms (forced to be the most specific for each gene) using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }else{
            message(sprintf("Second, pre-compute semantic similarity between %d terms using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=T)
        }
    }
    ## pre-compute semantic similarity between terms in subject
    sim.term <- suppressMessages(dDAGtermSim(ig, terms=terms, method=method.term, verbose=T))
    
    if(verbose){
        message(sprintf("Last, calculate pair-wise semantic similarity between %d genes using %s method (%s)...", length(genes), method.gene, as.character(Sys.time())), appendLF=T)
    }
    ## calculate pair-wise semantic similarity between input genes
    sim <- Matrix::Matrix(0, nrow=length(genes), ncol=length(genes), sparse=T)
    for(i in 1:(length(genes2terms)-1)){
        terms1 <- genes2terms[[i]]
        ind1 <- match(terms1, terms)
        
        progress_indicate(i, length(genes2terms), 10, flag=T)
        
        for(j in (i+1):length(genes2terms)){
            terms2 <- genes2terms[[j]]
            ind2 <- match(terms2, terms)
            
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            
            ## gene-gene similarity
            switch(method.gene, 
                average={res <- mean(sim12)},
                max={res <- max(sim12)},
                BM.average={Max_i2j <- apply(sim12,1,max); Max_j2i <- apply(sim12,2,max); res <- 0.5*(mean(Max_i2j) + mean(Max_j2i))},
                BM.max={Max_i2j <- apply(sim12,1,max); Max_j2i <- apply(sim12,2,max); res <- max(mean(Max_i2j), mean(Max_j2i))},
                BM.complete={Max_i2j <- apply(sim12,1,max); Max_j2i <- apply(sim12,2,max); res <- min(c(Max_i2j,Max_j2i))}
            )
            
            sim[i,j] <- res
            sim[j,i] <- res   
        }
    }
    rownames(sim) <- colnames(sim) <- genes
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    return(sim)
}

