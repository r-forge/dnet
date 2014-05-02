#' Function to estimate RWR dating based sample relationships from the input gene-sample matrix and graph
#'
#' \code{dRWRdating} is supposed to estmate sample relationships from the given gene-sample matrix and graph according to random walk restart (RWR) dating based method. The method includes: 1) RWR in the graph using the input matrix as seeds; 2) calculation of dating based contact strength (inner products of RWR-smooth columns of intput matrix); 3) estimation of the contact signficance by a randomalisation procedure.
#'
#' @param data an input gene-sample data matrix used for seeds
#' @param g an object of class "igraph" or "graphNEL"
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param num.permutation the number of permutations used to for generating the distribution of contact strength under randomalisation
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param adjp.cutoff the cutoff of adjusted pvalue to construct the contact graph
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' an object of class "dContact", a list with following components:
#' \itemize{
#'  \item{\code{zscore}: a symmetric matrix storing zscore between pairwise samples}
#'  \item{\code{pval}: a symmetric matrix storing pvalue between pairwise samples}
#'  \item{\code{adjpval}: a symmetric matrix storing adjusted pvalue between pairwise samples}

#'  \item{\code{cgraph}: the constructed contact graph (as a 'igraph' object) under the cutoff of adjusted value}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dRWR}}
#' @include dRWRdating.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#' V(subg)$name <- 1:vcount(subg)
#'
#' # 3) obtain the pre-computated affinity matrix
#' PTmatrix <- dRWR(subg, normalise="laplacian", restart=0.75)
#' # visualise affinity matrix
#' visHeatmapAdv(PTmatrix, Rowv=FALSE, Colv=FALSE, colormap="wyr", KeyValueName="Affinity")
#' 
#' # 3) obtain affinity matrix given sets of seeds
#' # define sets of seeds as data
#' # each seed with equal weight (i.e. all non-zero entries are '1')
#' aSeeds <- c(1,0,1,0,1)
#' bSeeds <- c(0,0,1,0,1)
#' data <- data.frame(aSeeds,bSeeds)
#' rownames(data) <- 1:5
#' # calcualte their two contact graph
#' dContact <- dRWRdating(data=data, g=subg)
#' dContact

dRWRdating <- function(data, g, normalise=c("laplacian","row","column","none"), restart=0.5, normalise.affinity.matrix=c("none","quantile"), num.permutation=10, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), adjp.cutoff=0.05, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    p.adjust.method <- match.arg(p.adjust.method)
    
    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
        if(ncol(data)<2){
            stop("The input data must be matrix with at least two columns.\n")
        }
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }
    if(is.null(rownames(data))) {
        stop("The function must require the row names of the input data.\n")
    }else if(any(is.na(rownames(data)))){
        warning("Data with NA as row names will be removed")
        data <- data[!is.na(rownames(data)),]
    }
    cnames <- colnames(data)
    if(is.null(cnames)){
        cnames <- seq(1,ncol(data))
    }
    ## check input graph
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
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
    
    if(verbose){
        message(sprintf("First, RWR on input graph (%d nodes and %d edges) using input matrix (%d rows and %d columns) as seeds (%s)...", vcount(ig), ecount(ig), nrow(data), ncol(data), as.character(Sys.time())), appendLF=T)
    }
    ## RWR of a Matrix
    PTmatrix <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, setSeeds=data, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix)))
    
    if(verbose){
        message(sprintf("Second, calculate dating-based contact strength (%s)...", as.character(Sys.time())), appendLF=T)
    }
    D <- t(PTmatrix)
    n <- nrow(D)
    obs <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
        obs[i, ] <- (D[i, ] %*% t(D))
    }
    rownames(obs) <- rownames(D)
    colnames(obs) <- rownames(D)
    
    B <- num.permutation
    if(verbose){
        message(sprintf("Third, generate the distribution of contact strength using %d permutations on nodes (%s)...", B, as.character(Sys.time())), appendLF=T)
    }
    exp_b <- list()
    for (b in 1:B){
        progress_indicate(b, B, 10, flag=T)
        seeds_random <- data[sample(1:nrow(data)),]
        rownames(seeds_random) <- rownames(data)
        PT_random <- suppressWarnings(suppressMessages(dRWR(g=ig, normalise=normalise, setSeeds=seeds_random, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix)))
        D <- t(PT_random)
        exp_random <- matrix(0, nrow=n, ncol=n)
        for (i in 1:n) {
            exp_random[i, ] <- (D[i, ] %*% t(D))
        }
        rownames(exp_random) <- rownames(D)
        colnames(exp_random) <- rownames(D)
        exp_b[[b]] <- exp_random
    }

    if(verbose){
        message(sprintf("Last, estimate the significance of contact strength: zscore, pvalue, and %s adjusted-pvalue (%s)...", p.adjust.method, as.character(Sys.time())), appendLF=T)
    }
    
    ## for zscore
    exp_mean <- matrix(0, ncol=n, nrow=n)
    exp_square <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        exp_mean <- exp_mean + exp_b[[b]]
        exp_square <- exp_square + exp_b[[b]]^2
    }
    exp_mean <- exp_mean/B
    exp_square <- exp_square/B
    exp_std <- sqrt(exp_square-exp_mean^2)
    zscore <- (obs-exp_mean)/exp_std
    
    ## for pvalue
    counts <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        counts <- counts + (obs < exp_b[[b]])
    }
    pval <- counts/B
    
    ## for adjusted pvalue
    adjpval <- pval
    ## lower part
    flag_lower <- lower.tri(pval, diag=T)
    adjpval[flag_lower] <- stats::p.adjust(pval[flag_lower], method=p.adjust.method[1])
    ## upper part
    flag_upper <- upper.tri(pval, diag=T)
    adjpval[flag_upper] <- stats::p.adjust(pval[flag_upper], method=p.adjust.method[1])

    if(verbose){
        message(sprintf("Also, construct the contact graph under the cutoff %1.1e of adjusted-pvalue (%s)...", adjp.cutoff, as.character(Sys.time())), appendLF=T)
    }
    flag <- adjpval < adjp.cutoff
    adjmatrix <- flag
    adjmatrix[flag] <- zscore[flag]
    cgraph <- igraph::graph.adjacency(adjmatrix, mode="undirected", weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)

    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    dContact <- list(zscore = zscore, 
                     pval = pval,
                     adjpval = adjpval, 
                     cgraph = cgraph, 
                     call = match.call(), 
                     method = "dnet")
    class(dContact) <- "dContact"
    
    invisible(dContact)
}
