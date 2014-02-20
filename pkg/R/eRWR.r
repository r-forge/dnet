#' Function to implement Random Walk with Restart (RWR) to pre-compute affinity matrix for the input graph
#'
#' \code{eRWR} is supposed to implement Random Walk with Restart (RWR) to pre-compute affinity matrix for for nodes in the input graph with respect to the starting node (loop over every node in the graph)
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for RWR
#' 
#' @return 
#' \itemize{
#'  \item{\code{PTmatrix}: affinity matrix with the dimension of n X n, where n is the number of nodes in the input graph. Columns stand for nodes starting to walk from, and rows for nodes ending up to. Therefore, a column for a starting node represents a steady-state affinity vector that the starting node will visit all the end nodes in the graph}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute assocaited
#' @export
#' @seealso \code{\link{eGraphInduce}}
#' @include eRWR.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- eGraphInduce(g, V(g), knn=0)
#'
#' # 3) calculate the affinity matrix
#' PTmatrix <- eRWR(subg, normalise="laplacian", restart=0.75)
#'
#' # 4) visualise affinity matrix
#' graphics::image(PTmatrix, col=visColormap("wyr")(64), zlim=c(0,1))

eRWR <- function(g, normalise=c("laplacian","row","column","none"), restart=0.75)
{
    
    normalise <- match.arg(normalise)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if ("weight" %in% list.vertex.attributes(ig)){
        adjM <- get.adjacency(ig, type="both", attr="weight", edges=F, names=T, sparse=F)
    }else{
        adjM <- get.adjacency(ig, type="both", attr=NULL, edges=F, names=T, sparse=F)
    }
    
    A <- adjM!=0
    if(normalise == "row"){
        D <- diag(apply(A,1,sum)^(-1))
        nadjM <- D %*% adjM
    }else if(normalise == "column"){
        D <- diag(apply(A,1,sum)^(-1))
        nadjM <- adjM %*% D
    }else if(normalise == "laplacian"){
        D <- diag(apply(A,1,sum)^(-0.5))
        nadjM <- D %*% adjM %*% D
    }else{
        nadjM <- adjM
    }
    
    #nig <- graph.adjacency(nadjM, mode=c("undirected"), weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
    ## update the vertex attributes
    nattr <- list.vertex.attributes(ig)
    for(attr in nattr){
        #nig <- set.vertex.attribute(nig, name=attr, index=V(nig), get.vertex.attribute(ig,attr))
    }
    ## update the edge attributes    
    eattr <- list.edge.attributes(ig)
    for(attr in eattr){
        if (!("weight" %in% attr)){
            #nig <- set.edge.attribute(nig, name=attr, index=E(nig), get.edge.attribute(ig,attr))
        }
    }
    
    ################## RWR
    ## restarting prob
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
        r <- 0.75
    }else if(restart>1 && restart<100){
        r <- restart/100
    }else{
        r <- restart
    }
    ## stopping critera
    stop_delta <- 1e-10   # L1 norm of successive estimates 'PT' below the threshold 'stop_delta'
    stop_step <- 100      # maximum steps of iterations
    
    P0matrix <- matrix(as.numeric(nadjM>0),nrow=nrow(nadjM),ncol=ncol(nadjM))
    PTmatrix <- matrix(0,nrow=nrow(P0matrix),   ncol=ncol(P0matrix))
    for(j in 1:ncol(P0matrix)){
        P0 <- as.matrix(P0matrix[,j],ncol=1)
        
        ## Initializing variables
        step <- 0
        delta <- 1
        
        PT <- P0
        ## Iterative update till convergence (delta<=1e-10)
        while (delta>stop_delta && step<=stop_step){
            PX <- (1-r) * nadjM %*% PT + r * P0
    
            # p-norm of v: sum((abs(v).p)^(1/p))
            delta <- sum(abs(PX-PT))
    
            PT <- PX
            step <- step+1
        }
        PTmatrix[,j] <- PT
    }
    
    invisible(PTmatrix)
}