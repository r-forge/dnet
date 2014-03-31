#' Function to generate a subgraph of a direct acyclic graph (DAG) induced by given vertices
#'
#' \code{dDagInduce} is supposed to produce a subgraph induced by given vertices, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph" or "graphNET" object, a list of the vertices of the graph, and the mode defining the paths to the root of DAG. The resultant subgraph inherits the class from the input one. The induced subgraph contains exactly the vertices of interest and their defined paths to the root of DAG. 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param nodes_query the vertices for which the calculation is performed
#' @param path.mode the mode of paths induced by nodes in query. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @return 
#' \itemize{
#'  \item{\code{subg}: an induced subgraph, an object of class "igraph" or "graphNEL"}
#' }
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{dDagInduce}}
#' @include dDagInduce.r
#' @examples
#' # 1) load GO Molelular Function as igraph object
#' load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
#' g <- ig.GOMF
#'
#' # 2) randomly select vertices as the query nodes
#' # the query nodes can be igraph vertex sequences
#' nodes_query <- V(g)[sample(V(g),5)]
#' # the more common, the query nodes can be term id
#' nodes_query <- V(g)[sample(V(g),5)]$name
#'
#' # 3) obtain the induced subgraph
#' # 3a) based on all possible paths (i.e. the complete subgraph induced)
#' subg <- dDagInduce(g, nodes_query, path.mode="all_paths")
#' # 3b) based on shortest paths (i.e. the most concise subgraph induced)
#' subg <- dDagInduce(g, nodes_query, path.mode="shortest_paths")

dDagInduce <- function (g, nodes_query, path.mode=c("all_paths","shortest_paths","all_shortest_paths"))
{
    
    path.mode <- match.arg(path.mode)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if(class(nodes_query)=="igraph.vs"){
        nodes_query <- nodes_query$name
    }
    
    ## check nodes in query
    ind <- match(nodes_query, V(ig)$name)
    nodes_mapped <- nodes_query[!is.na(ind)]
    if(length(nodes_mapped)==0){
        stop("Nodes in query cannot be found in the input graph.\n")
    }else{
        nodes_query <- V(ig)[nodes_mapped]
    }

    ## DAG being induced from nodes in query
    if(path.mode=="all_paths"){
        edgelist <- get.edgelist(ig, names=T)
            
        nodeLookUp <- new.env(hash=T, parent=emptyenv())
        ## A function to iterate to the root
        buildInducedGraph <- function(node) {
            ## For node already visited, nothing to do
            if (exists(node, envir=nodeLookUp, mode="logical", inherits=FALSE)){
                return(1)
            }
            ## put the node into nodeLookUp, and get parents
            assign(node, TRUE, envir=nodeLookUp)
            adjNodes <- edgelist[edgelist[,2]==node, 1]
            if (length(adjNodes)==0){
                return(2)
            }
            for (i in 1:length(adjNodes)){
                buildInducedGraph(adjNodes[i])
            }
            ## criteria for quite the loop/lookup
            return(0)
        }
        startNodes <- nodes_query$name
        lapply(startNodes, buildInducedGraph)
        nodeInduced <- ls(nodeLookUp)
            
    }else if(path.mode=="all_shortest_paths"){
        nodeRoot <- V(ig)[which(V(g)$term_distance==0)]
        aspaths <- get.all.shortest.paths(ig, from=nodeRoot, to=nodes_query)
        nodeInduced <- unique(unlist(aspaths$res))
    }else if(path.mode=="shortest_paths"){
        nodeRoot <- V(ig)[which(V(ig)$term_distance==0)]
        vpaths <- get.shortest.paths(ig, from=nodeRoot, to=nodes_query)
        nodeInduced <- unique(unlist(vpaths))
    }
    
    subg <- induced.subgraph(ig, vids=nodeInduced)
    
    if(class(g)=="graphNEL"){
        subg <- igraph.to.graphNEL(subg)
    }

    return(subg)
}