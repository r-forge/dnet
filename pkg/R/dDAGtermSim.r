#' Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dDAGtermSim} is supposed to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param terms the terms/nodes between which pair-wise semantic similarity is calculated. If NULL, all terms in the input dag will be used for calcluation, which is very prohibitively expensive!
#' @param method the method used to measure semantic similarity between input terms. It can be "Resnik" for information content (IC) of most informative information ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{http://webdocs.cs.ualberta.ca/~lindek/papers/sim.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - diference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param precomputed logical to indicate whether a sparse matrix of children x ancestors is pre-computed for the subsequent identification of common ancestors. By default, it sets to true. It is advisable to have enabled it when there are many input terms
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix containing pair-wise semantic similarity between input terms
#' @note When there are many input terms, 'precomputed' shoulb be true; otherwise keep it as FALSE for no more than a handful of input terms
#' @export
#' @import Matrix
#' @seealso \code{\link{dDAGinduce}}, \code{\link{dDAGancestor}}
#' @include dDAGtermSim.r
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
#' terms <- sample(V(dag)$name, 5)
#' sim <- dDAGtermSim(g=dag, terms=terms, method=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), precomputed=FALSE, verbose=TRUE)
#' sim
#' }

dDAGtermSim <- function (g, terms=NULL, method=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), precomputed=TRUE, verbose=TRUE)
{
    
    method <- match.arg(method)
    
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
    IC <- V(ig)$IC
    names(IC) <- V(ig)$name
    
    ## checking input terms
    terms <- terms[!is.na(terms)]
    if(is.null(terms) || is.na(terms)){
        terms <- V(ig)$name
    }else{
        flag <- terms %in% V(ig)$name
        if(sum(flag)!=0){
            terms <- terms[flag]
        }else{
            terms <- V(ig)$name
        }
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
    
    sim <- Matrix::Matrix(0, nrow=length(terms), ncol=length(terms), sparse=T)
    if(precomputed){
    
        if(verbose){
            message(sprintf("Calculate semantic similarity between %d terms using %s method with pre-computed option enabled (%s)...", length(terms), method, as.character(Sys.time())), appendLF=T)
        }
    
        ## pre-compute a sparse matrix of children x ancestors
        sCP <- dDAGancestor(ig, term1=NULL, term2=NULL, verbose=T)
        
        ind <- match(terms, rownames(sCP))
        
        for(i in 1:length(terms)){
            ind_i <- ind[i]
            
            progress_indicate(i, length(terms), 100, flag=T)
            
            for(j in i:length(terms)){
                ind_j <- ind[j]
            
                term1 <- terms[i]
                term2 <- terms[j]
                
                ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                tmp <- apply(sCP[c(ind_i,ind_j),],2,sum)
                ancestors <- names(which(tmp==2))
                mica <- ancestors[which.max(IC[ancestors])]
            
                if(method=="Resnik"){
                    res <- IC[mica]
                }else if(method=="Lin"){
                    res <- 2*IC[mica]/(IC[term1]+IC[term2])
                }else if(method=="Schlicker"){
                    res <- (2*IC[mica]/(IC[term1]+IC[term2])) * (1 - 10^(-IC[mica]))
                }else if(method=="Jiang"){
                    res <- 1- min(1, IC[term1]+IC[term2]-2*IC[mica])
                }else if(method=="Pesquita"){
                    ## graph information content similarity related to Tanimoto-Jacard index
                    ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                    allancestors <- names(which(tmp>=1))
                    res <- sum(IC[ancestors]) / sum(IC[allancestors])
                }
            
                sim[i,j] <- res
                sim[j,i] <- res
            
            }
        }
        
        
    }else{
        
        if(verbose){
            message(sprintf("Calculate semantic similarity between %d terms using %s method (%s)...", length(terms), method, as.character(Sys.time())), appendLF=T)
        }
        
        for(i in 1:length(terms)){
        
            progress_indicate(i, length(terms), 100, flag=T)
    
            for(j in i:length(terms)){
            
                term1 <- terms[i]
                term2 <- terms[j]
            
                ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                ancestors <- dDAGancestor(g=ig, term1, term2, verbose=F)
                mica <- ancestors[which.max(IC[ancestors])]
            
                if(method=="Resnik"){
                    res <- IC[mica]
                }else if(method=="Lin"){
                    res <- 2*IC[mica]/(IC[term1]+IC[term2])
                }else if(method=="Schlicker"){
                    res <- (2*IC[mica]/(IC[term1]+IC[term2])) * (1 - 10^(-IC[mica]))
                }else if(method=="Jiang"){
                    res <- 1- min(1, IC[term1]+IC[term2]-2*IC[mica])
                }else if(method=="Pesquita"){
                    ## graph information content similarity related to Tanimoto-Jacard index
                    ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                    subg <- dDAGinduce(g=ig, nodes_query=c(term1,term2), path.mode="all_paths")
                    allancestors <- V(subg)$name
                    res <- sum(IC[ancestors]) / sum(IC[allancestors])
                }
            
                sim[i,j] <- res
                sim[j,i] <- res
            
            }
        }
    }
    
    rownames(sim) <- colnames(sim) <- terms
    
    return(sim)
}

