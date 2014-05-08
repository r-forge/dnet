#' Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data
#'
#' \code{dDAGtermSim} is supposed to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param terms the terms/nodes between which pair-wise semantic similarity is calculated. If NULL, all terms in the input dag will be used for calcluation, which is very prohibitively expensive!
#' @param method the method used to measure semantic similarity between input terms. It can be "Resnik" for information content (IC) of most informative information ancestor (MICA) (see \url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{http://webdocs.cs.ualberta.ca/~lindek/papers/sim.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - diference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186})). By default, it uses "Schlicker" method
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts.
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns a sparse matrix containing pair-wise semantic similarity between input terms
#' @note none
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
#' sim <- dDAGtermSim(g=dag, terms=terms, method="Schlicker")
#' sim
#' }

dDAGtermSim <- function (g, terms=NULL, method=c("Resnik","Lin","Schlicker","Jiang","Pesquita")[3], fast=T, verbose=TRUE)
{
    
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
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
    
    if(verbose){
        message(sprintf("Calculate semantic similarity between %d terms using %s method (%s)...", length(terms), method, as.character(Sys.time())), appendLF=T)
    }
    
    ## pre-compute a sparse matrix of children x ancestors
    sCP <- dDAGancestor(ig, term1=terms, term2=NULL, verbose=T)
    
    allterms <- V(ig)$name
    ind <- match(terms, allterms)
    
    ## calculate pair-wise semantic similarity between input terms
    sim <- Matrix::Matrix(0, nrow=length(terms), ncol=length(terms), sparse=T)
    
    for(i in 1:length(terms)){
        ancestor_i <- which(sCP[i,]==1)
        
        progress_indicate(i, length(terms), 100, flag=T)
        
        if(fast){
            
            mat <- sCP[i:length(terms),]
            ancestor_js <- which(matrix(mat==1, nrow=length(terms)-i+1), arr.ind=T)
        
            flag <- is.element(ancestor_js[,2], ancestor_i)
            ca_js <- ancestor_js[flag,]
            ca_js_list <- split(ca_js[,2],ca_js[,1])
            
            mica_js <- sapply(ca_js_list, function(x){
                x[which.max(IC[x])]
            })
            js <- as.numeric(names(ca_js_list))+i-1
        
            if(method=="Resnik"){
                res <- IC[mica_js]
            }else if(method=="Lin"){
                res <- 2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])
            }else if(method=="Schlicker"){
                res <- (2*IC[mica_js]/(IC[ind[i]]+IC[ind[js]])) * (1 - 10^(-IC[mica_js]))
            }else if(method=="Jiang"){
                tmp <- IC[ind[i]]+IC[ind[js]]-2*IC[mica_js]
                tmp[tmp>1] <- 1
                res <- 1- tmp
            }else if(method=="Pesquita"){
                
                ## graph information content similarity related to Tanimoto-Jacard index
                ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                ## for all ancestors
                allan_js_list <- split(ancestor_js[,2],ancestor_js[,1])
                allan_union <- sapply(allan_js_list, function(x){
                    ux <- union(x, ancestor_i)
                    sum(IC[ux])
                })
                ## for all common ancestors
                allca_union <- sapply(ca_js_list, function(x){
                    sum(IC[x])
                })
                ## their ratio
                res <- allca_union / allan_union
                
            }
            
            sim[i,js] <- res
            sim[js,i] <- res
        
        }else{
            for(j in i:length(terms)){
                ancestor_j <- which(sCP[j,]==1)
                
                ## find the common ancestor with the maximal information content: most informative information ancestor (MICA)
                ancestors <- intersect(ancestor_i, ancestor_j)
                mica <- ancestors[which.max(IC[ancestors])]
            
                if(method=="Resnik"){
                    res <- IC[mica]
                }else if(method=="Lin"){
                    res <- 2*IC[mica]/(IC[ind[i]]+IC[ind[j]])
                }else if(method=="Schlicker"){
                    res <- (2*IC[mica]/(IC[ind[i]]+IC[ind[j]])) * (1 - 10^(-IC[mica]))
                }else if(method=="Jiang"){
                    res <- 1- min(1, IC[ind[i]]+IC[ind[j]]-2*IC[mica])
                }else if(method=="Pesquita"){
                    ## graph information content similarity related to Tanimoto-Jacard index
                    ## summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2
                    allancestors <- union(ancestor_i, ancestor_j)
                    res <- sum(IC[ancestors]) / sum(IC[allancestors])
                }
            
                sim[i,j] <- res
                sim[j,i] <- res
            }
        }
    }
    rownames(sim) <- colnames(sim) <- terms
    
    sim[is.na(sim)] <- 0
    
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
