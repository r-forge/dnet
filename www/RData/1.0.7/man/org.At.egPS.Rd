\name{org.At.egPS}
\alias{org.At.egPS}

\title{Annotations of Arabidopsis Entrez Genes (EG) by phylostratific age (PS).}


\usage{
org.At.egPS <- dRDataLoader(RData='org.At.egPS')
}

\description{
An R object that contains phylostratific age information for Arabidopsis Entrez Genes. This data is prepared based on 1) SUPERFAMILY database which provides domain architecture assignments to all completely sequenced genomes including eukaryotic genomes; 2) ancestral domain architecture repertoires inferred by applying Dollo parsimony to eukaryotic part of species tree of life (sTOL), from which the most recent common ancestor of each domain architecture is determined. The domain architecture for an Entrez gene is the protein product with the longest length of amino acids. Thus, phylostratific age for a Arabidopsis Entrez gene is the first appearance of its domain architecture along the branch from the eukaryotic ancestor to the arabidopsis, and thus can be measured by: i) the most recent common ancestor, ii) how many steps it is away starting from the eukaryotic ancestor, and how far it is in the terms of the branch length from the eukaryotic ancestor. 
}

\value{
    an object of class "GS", a list with following components:
    \itemize{
        \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene sets (i.e. phylogenetic placement along the branch starting from the eukaryotic ancestor). The 4 columns are "setID" (i.e. "phylogenetic placement ID"), "name" (i.e. name for that placement in the form of "TaxonID:Name"), "namespace" (i.e. Rank for that placement) and "distance" (i.e. the branch length from the eukaryotic ancestor). Notably, since the sTOL is bifurcating with exactly two descendants (unlike the multifurcating nature of the NCBI taxonomy), an internal node in sTOL is either mapped onto a unique taxonomic identifier or left empty (assumedly a hypothetical unknown ancestor). In the latter case, hypothetical unknown ancestor is filled with the information in its nearest descendant with known taxonomic information.
        }
        \item{\code{gs}: a list of gene sets, each storing gene members thereof. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
    }
}

\references{
Morais et al. (2011) SUPERFAMILY 1.75 including a domain-centric gene ontology method. \emph{Nucleic Acids Res}, 39(Database issue):D427-34.
\cr
Fang et al. (2013) A daily-updated tree of (sequenced) life as a reference for genome research. \emph{Scientific reports}, 3:2015.
}

\keyword{datasets}

\examples{
org.At.egPS <- dRDataLoader(RData='org.At.egPS')
names(org.At.egPS)
}
