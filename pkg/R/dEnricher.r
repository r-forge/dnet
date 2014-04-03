#' Function to conduct enrichment analysis given the input data and the ontology in query
#'
#' \code{dEnricher} is supposed to conduct enrichment analysis given the input data and the ontology in query. It returns an object of class "eTerm". 
#'
#' @param data an input vector. It contains either Entrez Gene ID or Symbol
#' @param identity the type of gene identity (i.e. row names of input data), either "symbol" for gene symbols (by default) or "entrez" for Entrez Gene ID. The option "symbol" is preferred as it is relatively stable from one update to another; also it is possible to search against synonyms (see the next parameter)
#' @param check.symbol.identity logical to indicate whether synonyms will be searched against when gene symbols cannot be matched. By default, it sets to FALSE since it may take a while to do such check using all possible synoyms
#' @param genome the genome identity. It can be one of "Hs" for human, "Mm" for mouse, "Rn" for rat, "Gg" for chicken, "Ce" for c.elegans, "Dm" for fruitfly, "Da" for zebrafish, and "At" for arabidopsis
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, and the molecular signatures database (Msigdb) in human (including "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7"). Note: These four ("GOBP", "GOMF", "GOCC" and "PS") are availble for all genomes/species; for "Hs" and "Mm", these five ("DO", "HPPA", "HPMI", "HPON" and "MP") are also supported; all "Msigdb" are only supported in "Hs". For details on the eligibility for pairs of input genome and ontology, please refer to the online Documentations at \url{http://dnet.r-forge.r-project.org/docs.html}
#' @param sizeRange the minimum and maximum size of members of each gene set in consideration. By default, it sets to a minimum of 10 but no more than 1000
#' @param which_distance which distance of terms in the ontology is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the statistic test used. It can be "FisherTest" for using fisher's exact test, or "HypergeoTest" for using hypergeometric test. Hypergeometric test is to sample at random from the background containing annotated and non-annotated genes, and thus compare sampling to background. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to respect the hierarchy of the ontology. It can be "none" for not considering the ontology heirarchy, "elim" for the alogrithm to computer the significance of a term in terms of the significance of its children (precisely, once genes are already annotated to a signficantly enriched term, all these genes are eliminated from the ancestors of that term)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at \url{"http://dnet.r-forge.r-project.org/data"}. For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. Surely, the location can be anywhere as long as the user provides the correct path pointing to (otherwise, the script will have to remote download each time). Here is the UNIX command for downloading all RData files (preserving the directory structure): \eqn{wget -r -l2 -A "*.RData" -np -nH --cut-dirs=0 "http://dnet.r-forge.r-project.org/data"}
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene set in consideration, and the 4 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{gs}: a list of gene sets, each storing gene members. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note None
#' @export
#' @seealso \code{\link{dEnricher}}
#' @include dEnricher.r
#' @examples
#' \dontrun{
#' load(url("http://dnet.r-forge.r-project.org/data/Datasets/Hiratani_TableS1.RData"))
#' data <- rownames(RT)[1:1000]
#' eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="MP", RData.location="./RData_Rd")
#' eTerm <- dEnricher(data, identity="symbol", genome="Mm", ontology="MP", ontology.algorithm="elim", RData.location="./RData_Rd")
#' cbind(eTerm$set_info[which(eTerm$pvalue < 1e-3), c(1,2)], eTerm$pvalue[which(eTerm$pvalue < 1e-3)])
#'
#' # highlight the top significant terms and also color-code all terms according to the adjust p-values
#' load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.MP.RData"))
#' g <- ig.MP
#' nodes_query <- names(sort(eTerm$adjp)[1:5])
#' nodes.highlight <- rep("red", length(nodes_query))
#' names(nodes.highlight) <- nodes_query
#' subg <- dDAGinduce(g, nodes_query)
#' visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,2), node.attrs=list(color=nodes.highlight))
#' }

dEnricher <- function(data, identity=c("symbol","entrez"), check.symbol.identity=FALSE, genome=c("Hs", "Mm", "Rn", "Gg", "Ce", "Dm", "Da", "At"), ontology=c("GOBP","GOMF","GOCC","PS","DO","HPPA","HPMI","HPON","MP", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7"), sizeRange=c(10,1000), which_distance=NULL, test=c("FisherTest","HypergeoTest"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","elim"), verbose=T, RData.location="http://dnet.r-forge.r-project.org/data")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    identity <- match.arg(identity)
    genome <- match.arg(genome)
    ontology <- match.arg(ontology)
    test <- match.arg(test)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    
    if (is.vector(data)){
        data <- unique(data)
    }else{
        stop("The input data must be a vector.\n")
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("First, load the ontology %s and its gene associations in the genome %s (%s) ...", ontology, genome, as.character(now)), appendLF=T)
    }
    
    ###############################
    ## check the eligibility for pairs of input genome and ontology
    all.ontologies <- c("GOBP","GOMF","GOCC","PS","DO","HPPA","HPMI","HPON","MP", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
    possible.ontologies <- switch(genome,
                       Hs = all.ontologies[c(1:4, 5:9, 10:24)],
                       Mm = all.ontologies[c(1:4, 5:9)],
                       Rn = all.ontologies[c(1:4)],
                       Gg = all.ontologies[c(1:4)],
                       Ce = all.ontologies[c(1:4)],
                       Dm = all.ontologies[c(1:4)],
                       Da = all.ontologies[c(1:4)],
                       At = all.ontologies[c(1:4)]
                       )
    if(!(ontology %in% possible.ontologies)){
        stop(sprintf("The input pair of genome (%s) and ontology (%s) are not supported.\nThe supported ontologies in genome (%s): %s.\n", genome, ontology, genome, paste(possible.ontologies,collapse=", ")))
    }

    ###############################
    ## make sure there is no "/" at the end
    path_host <- gsub("/$", "", RData.location)
    if(path_host=="" || length(path_host)==0 || is.na(path_host)){
        path_host <- "http://dnet.r-forge.r-project.org/data"
    }

    #########
    ## load Enterz Gene information
    EG <- list()
    load_EG_remote <- paste(path_host, "/", genome, "/org.", genome, ".eg.RData", sep="")
    load_EG_local1 <- file.path(path_host, paste("data/", genome, "/org.", genome, ".eg.RData", sep=""))
    load_EG_local2 <- file.path(path_host, paste(genome, "/org.", genome, ".eg.RData", sep=""))
    load_EG_local3 <- file.path(path_host, paste("org.", genome, "/org.", genome, ".eg.RData", sep=""))
    ## first, load local R files
    EG_local <- c(load_EG_local1, load_EG_local2, load_EG_local3)
    load_flag <- sapply(EG_local, function(x){
        if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
        ifelse(file.exists(x), TRUE, FALSE)
    })
    ## otherwise, load remote R files
    if(sum(load_flag)==0){
        if(class(try(load(url(load_EG_remote)), T))=="try-error"){
            load_EG_remote <- paste("http://dnet.r-forge.r-project.org/data/", genome, "/org.", genome, ".eg.RData", sep="")
            if(class(try(load(url(load_EG_remote)), T))=="try-error"){
                stop("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")
            }
        }
        load_EG <- load_EG_remote
    }else{
        load_EG <- EG_local[load_flag]
        load(load_EG)
    }
    eval(parse(text=paste("EG <- org.", genome, ".eg", sep="")))
    
    if(verbose){
        message(sprintf("\tLoad Enterz Gene information from %s", load_EG), appendLF=T)
    }
    
    #########
    ## load annotation information
    GS <- list()
    load_GS_remote <- paste(path_host, "/", genome, "/org.", genome, ".eg", ontology, ".RData", sep="")
    load_GS_local1 <- file.path(path_host, paste("data/", genome, "/org.", genome, ".eg", ontology, ".RData", sep=""))
    load_GS_local2 <- file.path(path_host, paste(genome, "/org.", genome, ".eg", ontology, ".RData", sep=""))
    load_GS_local3 <- file.path(path_host, paste("org.", genome, ".eg", ontology, ".RData", sep=""))
    ## first, load local R files
    GS_local <- c(load_GS_local1, load_GS_local2, load_GS_local3)
    load_flag <- sapply(GS_local, function(x){
        if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
        ifelse(file.exists(x), TRUE, FALSE)
    })
    ## otherwise, load remote R files
    if(sum(load_flag)==0){
        if(class(try(load(url(load_GS_remote)), T))=="try-error"){
            load_GS_remote <- paste("http://dnet.r-forge.r-project.org/data/", genome, "/org.", genome, ".eg", ontology, ".RData", sep="")
            if(class(try(load(url(load_GS_remote)), T))=="try-error"){
                stop("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")
            }
        }
        load_GS <- load_GS_remote
    }else{
        load_GS <- GS_local[load_flag]
        load(load_GS)
    }
    eval(parse(text=paste("GS <- org.", genome, ".eg", ontology, sep="")))
    
    if(verbose){
        message(sprintf("\tLoad annotation information from %s", load_GS), appendLF=T)
    }
    
    ###############################
    allGeneID <- EG$gene_info$GeneID
    allSymbol <- as.vector(EG$gene_info$Symbol)
    allSynonyms <- as.vector(EG$gene_info$Synonyms)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Then, do mapping based on %s (%s) ...", identity, as.character(now)), appendLF=T)
    }
    
    if(identity == "symbol"){
    
        Symbol <- data
        
        ## correct for those symbols being shown as DATE format
        if(1){
            ## for those starting with 'Mar' in a excel-input date format
            a <- Symbol
            flag <- grep("-Mar$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Mar$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("March",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }

            ## for those starting with 'Sep' in a excel-input date format
            a <- Symbol
            flag <- grep("-Sep$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Sep$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("Sept",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }
        }
        
        ## case-insensitive
        match_flag <- match(tolower(Symbol),tolower(allSymbol))
        
        ## match vis Synonyms for those unmatchable by official gene symbols
        if(check.symbol.identity){
            ## match Synonyms (if not found via Symbol)
            na_flag <- is.na(match_flag)
            a <- Symbol[na_flag]

            ###
            tmp_flag <- is.na(match(tolower(allSymbol), tolower(Symbol)))
            tmp_Synonyms <- allSynonyms[tmp_flag]
            Orig.index <- seq(1,length(allSynonyms))
            Orig.index <- Orig.index[tmp_flag]
            ###

            b <- sapply(1:length(a), function(x){
                tmp_pattern1 <- paste("^",a[x],"\\|", sep="")
                tmp_pattern2 <- paste("\\|",a[x],"\\|", sep="")
                tmp_pattern3 <- paste("\\|",a[x],"$", sep="")
                tmp_pattern <- paste(tmp_pattern1,"|",tmp_pattern2,"|",tmp_pattern3, sep="")
                tmp_result <- grep(tmp_pattern, tmp_Synonyms, ignore.case=T, perl=T, value=F)
                ifelse(length(tmp_result)==1, Orig.index[tmp_result[1]], NA)
            })
            match_flag[na_flag] <- b
            
            if(verbose){
                now <- Sys.time()
                message(sprintf("\tAmong %d symbols of input data, there are %d mappable via official gene symbols, %d mappable via gene alias, but %d left unmappable", length(Symbol), (length(Symbol)-length(a)), sum(!is.na(b)), sum(is.na(b))), appendLF=T)
            }
        }else{
            if(verbose){
                now <- Sys.time()
                message(sprintf("\tAmong %d symbols of input data, there are %d mappable via official gene symbols but %d left unmappable", length(Symbol), (sum(!is.na(match_flag))), (sum(is.na(match_flag)))), appendLF=T)
            }
        
        }
        
        ## convert into GeneID
        GeneID <- allGeneID[match_flag]
        
    }else{
        GeneID <- data
        match_flag <- match(GeneID,allGeneID)
        GeneID <- allGeneID[match_flag]
    }
    
    genes.group <- GeneID[!is.na(GeneID)]
    
    ## filter based on "which_distance"
    if(!is.null(which_distance) & sum(is.na(GS$set_info$distance))==0){
        set_filtered <- sapply(which_distance, function(x) {
            GS$set_info$setID[(GS$set_info$distance==as.integer(x))]
        })
        set_filtered <- unlist(set_filtered)
    }else{
        set_filtered <- names(GS$gs)
    }
    ind.distance <- match(set_filtered,names(GS$gs))
    
    ## derive the "gs" of interest
    gs.length <- sapply(GS$gs, length)
    ind.length <- which(gs.length >= sizeRange[1] & gs.length <= sizeRange[2])
    
    ind <- intersect(ind.distance, ind.length)
    gs <- GS$gs[ind]
    set_info <- GS$set_info[ind,]
    nSet <- length(gs)
    
    if(nSet==0){
        stop("There is no gene set being used.\n")
    }

    ##############################################################################################
    # the statistical tests used to compare a set of genes of interest to a set of reference genes
    # two distinct ways to model the problem:
    # Hypergeometric test: sampling at random from the background containing annotated and non-annotated genes (hypergeometric test); thus compare sampling to background
    # Fisher's exact test: testing the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not); thus compare sampling to the left part of background (after sampling without replacement)
    doFisherTest <- function(genes.group, genes.term, genes.universe){
        genes.hit <- intersect(genes.group, genes.term)
        # num of success in sampling
        X <- length(genes.hit)
        # num of sampling
        K <- length(genes.group)
        # num of success in background
        M <- length(genes.term)
        # num in background
        N <- length(genes.universe)
        ## Prepare a two-dimensional contingency table: #success in sampling, #success in background, #failure in sampling, and #failure in left part
        cTab <- matrix(c(X, K-X, M-X, N-M-K+X), nrow=2, dimnames=list(c("anno", "notAnno"), c("group", "notGroup")))
        p.value <- ifelse(all(cTab==0), 1, stats::fisher.test(cTab, alternative="greater")$p.value)
        return(p.value)
    }

    doHypergeoTest <- function(genes.group, genes.term, genes.universe){
        genes.hit <- intersect(genes.group, genes.term)
        # num of success in sampling
        X <- length(genes.hit)
        # num of sampling
        K <- length(genes.group)
        # num of success in background
        M <- length(genes.term)
        # num in background
        N <- length(genes.universe)
    
        x <- X
        m <- M
        n <- N-M # num of failure in background
        k <- K
        p.value <- ifelse(m==0 | k==0, 1, stats::phyper(x,m,n,k, lower.tail=F, log.p=F))
        return(p.value)
    }
    ##############################################################################################
    
    ## force use classic ontology.algorithm when the ontology is derived from "Msigdb" or "PS"
    if(length(grep("Msigdb",ontology))>0 || ontology=="PS"){
        ontology.algorithm <- "classic"
    }
    
    terms <- names(gs)
    genes.universe <- as.numeric(unique(unlist(gs[terms])))
    genes.group <- intersect(genes.universe, genes.group)
    
    if(length(genes.group)==0){
        stop("There is no gene being used.\n")
    }
    
    if(ontology.algorithm=="none"){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("Third, perform enrichment analysis using %s (%s) ...", test, as.character(now)), appendLF=T)
            if(is.null(which_distance)){
                message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations", length(terms), paste(sizeRange,collapse=",")), appendLF=T)
            }else{
                message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations and [%s] distance", length(terms), paste(sizeRange,collapse=","), paste(which_distance,collapse=",")), appendLF=T)
            }
        }
    
        pvals <- sapply(terms, function(term){
            genes.term <- as.numeric(unique(unlist(gs[term])))
            p.value <- switch(test,
                FisherTest =  doFisherTest(genes.group, genes.term, genes.universe),
                HypergeoTest = doHypergeoTest(genes.group, genes.term, genes.universe)
            )
        })

    }else if(ontology.algorithm=="elim"){

        if(verbose){
            now <- Sys.time()
            message(sprintf("Third, perform enrichment analysis using %s based on %s algorithm to respect ontology structure (%s) ...", test, ontology.algorithm, as.character(now)), appendLF=T)
        }
    
        ###############################
        ## load ontology information
        g <- list()
        load_g_remote <- paste(path_host, "/Obo/ig.", ontology, ".RData", sep="")
        load_g_local1 <- file.path(path_host, paste("data/Obo/ig.", ontology, ".RData", sep=""))
        load_g_local2 <- file.path(path_host, paste("Obo/ig.", ontology, ".RData", sep=""))
        load_g_local3 <- file.path(path_host, paste("ig.", ontology, ".RData", sep=""))
        ## first, load local R files
        g_local <- c(load_g_local1, load_g_local2, load_g_local3)
        load_flag <- sapply(g_local, function(x){
            if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
            ifelse(file.exists(x), TRUE, FALSE)
        })
        ## otherwise, load remote R files
        if(sum(load_flag)==0){
            if(class(try(load(url(load_g_remote)), T))=="try-error"){
                load_g_remote <- paste("http://dnet.r-forge.r-project.org/data/Obo/ig.", ontology, ".RData", sep="")
                if(class(try(load(url(load_g_remote)), T))=="try-error"){
                    stop("Built-in Rdata files cannot be loaded. Please check your internet connection or Rdata location in your local machine.\n")
                }
            }
            load_g <- load_g_remote
        }else{
            load_g <- g_local[load_flag]
            load(load_g)
        }
        eval(parse(text=paste("g <- ig.", ontology, sep="")))
        
        if(verbose){
            message(sprintf("\tLoad ontology information from %s", load_g), appendLF=T)
        }
        
        ###############################
        subg <- dDAGinduce(g, terms, path.mode="all_paths")
        
        if(verbose){
            message(sprintf("\tThere are %d terms being used", length(V(subg))), appendLF=T)
        }
        
        level2node <- dDAGlevel(subg, level.mode="longest_path", return.mode="level2node")
        
        ## build a hash environment from the named list "level2node"
        ## level2node.Hash: key (level), value (a list of nodes/terms)
        level2node.Hash <- list2env(level2node)
        ## ls(level2node.Hash)

        ## create a new (empty) hash environment
        ## sigNode2pval.Hash: key (node called significant), value (pvalue)
        sigNode2pval.Hash <- new.env(hash=T, parent=emptyenv())
        ## node2pval.Hash: key (node at ancestor), value (genes to be eliminated)
        ancNode2gene.Hash <- new.env(hash=T, parent=emptyenv())
        ## node2pval.Hash: key (node), value (pvalue)
        node2pval.Hash <- new.env(hash=T, parent=emptyenv())

        nLevels <- length(level2node)
        pval.cutoff <- 0.01
        for(i in nLevels:1) {
            currNodes <- get(as.character(i), envir=level2node.Hash, mode="character")
            currAnno <- GS$gs[currNodes]
    
            ## update "ancNode2gene.Hash" for each node/term
            for(currNode in currNodes){
                genes.term <- unique(unlist(GS$gs[currNode]))
        
                ## remove the genes (if any already marked) from annotations by the current node/term
                if(exists(currNode, envir=ancNode2gene.Hash, mode="numeric")){
                    genes.elim <- get(currNode, envir=ancNode2gene.Hash, mode="numeric")
                    genes.term <- setdiff(genes.term, genes.elim)
                    #message(sprintf("\t\t%d %d", length(genes.elim), length(genes.term)), appendLF=T)
                }
        
                ## do fisher exact test
                pvalue <- doFisherTest(genes.group, genes.term, genes.universe)
        
                ## store the result (the p-value)
                assign(currNode, pvalue, envir=node2pval.Hash)
        
                ## condition to update "ancNode2gene.Hash"
                if(pvalue < pval.cutoff) {
                    ## mark the significant node
                    assign(currNode, pvalue, envir=sigNode2pval.Hash)

                    ## retrieve genes annotated by the significant node for the subsequent eliminating
                    elimGenesID <- currAnno[[currNode]]

                    ## find all the ancestors of the significant node
                    dag.ancestors <- dDAGinduce(subg, currNode, path.mode="all_paths")
                    ancestors <- setdiff(V(dag.ancestors)$name, currNode)
            
                    ## get only those ancestors that are already in "ancNode2gene.Hash"
                    oldAncestors2GenesID <- sapply(ancestors, function(ancestor){
                        if(exists(ancestor, envir=ancNode2gene.Hash, mode="numeric")){
                            get(ancestor, envir=ancNode2gene.Hash, mode='numeric')
                        }
                    })

                    ## add the new GenesID to the ancestors
                    newAncestors2GenesID <- lapply(oldAncestors2GenesID, function(oldGenes){
                        union(oldGenes, elimGenesID)
                    })

                    ## update the "ancNode2gene.Hash" table
                    if(length(newAncestors2GenesID) > 0){
                        sapply(names(newAncestors2GenesID), function(ancestor){
                            assign(ancestor, newAncestors2GenesID[[ancestor]], envir=ancNode2gene.Hash)
                        })
                    }
                }
            }
    
            if(verbose){
                num.signodes <- length(ls(sigNode2pval.Hash))
                num.ancnodes <- length(ls(ancNode2gene.Hash))
                num.elimgenes <- length(unique(unlist(as.list(ancNode2gene.Hash))))
                message(sprintf("\tAt level %d, there are %d nodes/terms: so far, a total of %d nodes called significant, and %d ancestral nodes changed (%d genes eliminated)", i, length(currNodes), num.signodes, num.ancnodes, num.elimgenes), appendLF=T)
            }
        }
        pvals <- unlist(as.list(node2pval.Hash))
    
    }

    if(verbose){
        now <- Sys.time()
        message(sprintf("Last, adjust the p-values using the %s method (%s) ...", p.adjust.method, as.character(now)), appendLF=T)
    }

    ## Adjust P-values for multiple comparisons
    adjpvals <- stats::p.adjust(pvals, method=p.adjust.method)
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    eTerm <- list(set_info = set_info,
                  gs       = gs,
                  data     = data,
                  pvalue   = pvals,
                  adjp     = adjpvals,
                  call     = match.call()
                 )
    class(eTerm) <- "eTerm"
    
    invisible(eTerm)
}