#' Function to load dnet built-in RData
#'
#' \code{dRDataLoader} is supposed to load dnet built-in RData.
#'
#' @param RData which built-in RData to load. It can be one of "TCGA_mutations", "ig.DO", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.At.eg", "org.At.egGOBP", "org.At.egGOCC", "org.At.egGOMF", "org.At.egPS", "org.At.egSF", "org.At.string", "org.Ce.eg", "org.Ce.egGOBP", "org.Ce.egGOCC", "org.Ce.egGOMF", "org.Ce.egPS", "org.Ce.egSF", "org.Ce.string", "org.Da.eg", "org.Da.egGOBP", "org.Da.egGOCC", "org.Da.egGOMF", "org.Da.egPS", "org.Da.egSF", "org.Da.string", "org.Dm.eg", "org.Dm.egGOBP", "org.Dm.egGOCC", "org.Dm.egGOMF", "org.Dm.egPS", "org.Dm.egSF", "org.Dm.string", "org.Gg.eg", "org.Gg.egGOBP", "org.Gg.egGOCC", "org.Gg.egGOMF", "org.Gg.egPS", "org.Gg.egSF", "org.Gg.string", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.string", "org.Mm.eg", "org.Mm.egDO", "org.Mm.egGOBP", "org.Mm.egGOCC", "org.Mm.egGOMF", "org.Mm.egHPCM", "org.Mm.egHPMA", "org.Mm.egHPMI", "org.Mm.egHPPA", "org.Mm.egMP", "org.Mm.egPS", "org.Mm.egSF", "org.Mm.string", "org.Rn.eg", "org.Rn.egGOBP", "org.Rn.egGOCC", "org.Rn.egGOMF", "org.Rn.egPS", "org.Rn.egSF", "Hiratani_TableS1", "CLL", "org.Rn.string". On the meanings, please refer to the Documentations
#' @param genome the genome identity. It can be one of "Hs" for human, "Mm" for mouse, "Rn" for rat, "Gg" for chicken, "Ce" for c.elegans, "Dm" for fruitfly, "Da" for zebrafish, and "At" for arabidopsis
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for domain superfamily assignments, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, and Drug-Gene Interaction database (DGIdb) and the molecular signatures database (Msigdb) only in human (including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7"). Note: These four ("GOBP", "GOMF", "GOCC" and "PS") are availble for all genomes/species; for "Hs" and "Mm", these six ("DO", "HPPA", "HPMI", "HPCM", "HPMA" and "MP") are also supported; all "Msigdb" are only supported in "Hs". For details on the eligibility for pairs of input genome and ontology, please refer to the online Documentations at \url{http://supfam.org/dnet/docs.html}
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at \url{http://dnet.r-forge.r-project.org/data}. For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. Surely, the location can be anywhere as long as the user provides the correct path pointing to (otherwise, the script will have to remotely download each time). Here is the UNIX command for downloading all RData files (preserving the directory structure): \eqn{wget -r -l2 -A "*.RData" -np -nH --cut-dirs=0 "http://dnet.r-forge.r-project.org/data"}
#' @return 
#' any use-specified variable that is given on the right side of the assigement sign '<-', which contains the loaded RData.
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no RData will be loaded onto the working environment. 
#' @export
#' @seealso \code{\link{dRDataLoader}}
#' @include dRDataLoader.r
#' @examples
#' \dontrun{
#' org.Hs.egSF <- dRDataLoader(RData='org.Hs.egSF')
#' org.Hs.eg <- dRDataLoader(RData='org.Hs.eg')
#' org.Hs.egDGIdb <- dRDataLoader(RData='org.Hs.egDGIdb')
#' org.Hs.egMsigdbC2KEGG <- dRDataLoader(RData='org.Hs.egMsigdbC2KEGG')
#' org.Hs.egHPPA <- dRDataLoader(genome='Hs', ontology='HPPA')
#' ig.MP <- dRDataLoader(RData='ig.MP')
#' }

dRDataLoader <- function(RData=c(NA,"TCGA_mutations", "ig.DO", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.At.eg", "org.At.egGOBP", "org.At.egGOCC", "org.At.egGOMF", "org.At.egPS", "org.At.egSF", "org.At.string", "org.Ce.eg", "org.Ce.egGOBP", "org.Ce.egGOCC", "org.Ce.egGOMF", "org.Ce.egPS", "org.Ce.egSF", "org.Ce.string", "org.Da.eg", "org.Da.egGOBP", "org.Da.egGOCC", "org.Da.egGOMF", "org.Da.egPS", "org.Da.egSF", "org.Da.string", "org.Dm.eg", "org.Dm.egGOBP", "org.Dm.egGOCC", "org.Dm.egGOMF", "org.Dm.egPS", "org.Dm.egSF", "org.Dm.string", "org.Gg.eg", "org.Gg.egGOBP", "org.Gg.egGOCC", "org.Gg.egGOMF", "org.Gg.egPS", "org.Gg.egSF", "org.Gg.string", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.string", "org.Mm.eg", "org.Mm.egDO", "org.Mm.egGOBP", "org.Mm.egGOCC", "org.Mm.egGOMF", "org.Mm.egHPCM", "org.Mm.egHPMA", "org.Mm.egHPMI", "org.Mm.egHPPA", "org.Mm.egMP", "org.Mm.egPS", "org.Mm.egSF", "org.Mm.string", "org.Rn.eg", "org.Rn.egGOBP", "org.Rn.egGOCC", "org.Rn.egGOMF", "org.Rn.egPS", "org.Rn.egSF", "Hiratani_TableS1", "CLL", "org.Rn.string"), genome=c(NA, "Hs", "Mm", "Rn", "Gg", "Ce", "Dm", "Da", "At"), ontology=c(NA, "GOBP","GOMF","GOCC","PS","PS2","SF","DO","HPPA","HPMI","HPCM","HPMA","MP", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb"), RData.location="http://supfam.org/dnet/data")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    RData <- match.arg(RData)
    genome <- match.arg(genome)
    ontology <- match.arg(ontology)
    
    ###############################
    if(is.na(RData)){
        if(is.na(genome) | is.na(ontology)){
            stop(sprintf("Please make sure that either 'RData' is NOT NA or both 'genome' and 'ontology' are NOT NA.\n"))
        }else{
            
            ###############################
            ## get RData for annotations from input 'genome' and 'ontology'
            ###############################
                               
			## check the eligibility for pairs of input genome and ontology
			all.ontologies <- c("GOBP","GOMF","GOCC","PS","PS2","SF","DO","HPPA","HPMI","HPCM","HPMA","MP", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb")
			possible.ontologies <- switch(genome,
							   Hs = all.ontologies[c(1:6, 7:12, 13:28, 29)],
							   Mm = all.ontologies[c(1:6, 7:12)],
							   Rn = all.ontologies[c(1:6)],
							   Gg = all.ontologies[c(1:6)],
							   Ce = all.ontologies[c(1:6)],
							   Dm = all.ontologies[c(1:6)],
							   Da = all.ontologies[c(1:6)],
							   At = all.ontologies[c(1:6)]
							   )
			if(!(ontology %in% possible.ontologies)){
				stop(sprintf("The input pair of genome (%s) and ontology (%s) are not supported.\nThe supported ontologies in genome (%s): %s.\n", genome, ontology, genome, paste(possible.ontologies,collapse=", ")))
            }else{
                RData <- paste(paste('org.', genome, '.eg', ontology, sep=''))
            }
        }
    }
    
    ###############################
    
    ###############################
    ## make sure there is no "/" at the end
    path_host <- gsub("/$", "", RData.location)
    if(path_host=="" || length(path_host)==0 || is.na(path_host)){
        path_host <- "http://dnet.r-forge.r-project.org/data"
    }
    
    ## load 
    load_remote <- paste(path_host, "/", RData, ".RData", sep="")
    load_local1 <- file.path(path_host, paste("data/", RData, ".RData", sep=""))
    load_local2 <- file.path(path_host, paste(RData, ".RData", sep=""))
    load_package <- RData
    
    
    ## first, load data from the package itself
    if(length(suppressWarnings(tryCatch(eval(parse(text=paste("data(",load_package,", package='dnet')",sep=""))), error=function(e) e, warning=function(w) w)))==2){
        ## second, load local R files
        RData_local <- c(load_local1, load_local2)
        load_flag <- sapply(RData_local, function(x){
            if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
            ifelse(file.exists(x), TRUE, FALSE)
        })
        ## otherwise, load remote R files
        if(sum(load_flag)==0){
        
            con <- url(load_remote)
            if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                load_remote <- paste("http://supfam.org/dnet/data/", RData, ".RData", sep="")
            
                con <- url(load_remote)
                if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                    load_remote <- paste("http://dnet.r-forge.r-project.org/data/", RData, ".RData", sep="")
                
                    con <- url(load_remote)
                    if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                        load_remote <- paste("https://github.com/hfang-bristol/dnet/tree/master/data/", RData, ".RData", sep="")
                        
                        con <- url(load_remote)
                        if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                            stop("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")   
                        }
                    }
                }
            }
            close(con)
        
            load_RData <- load_remote
        }else{
            load_RData <- RData_local[load_flag]
            load(load_RData)
        }
    }else{
        load_RData <- sprintf("package 'dnet' version %s", packageVersion("dnet"))
    }
    
    out <- ''
    eval(parse(text=paste("out <- ", RData, sep="")))
    message(sprintf("'%s' (from %s) has been loaded into the working environment", RData, load_RData), appendLF=T)
    
    invisible(out)
}