\name{org.Gg.egSF}
\alias{org.Gg.egSF}

\title{Annotations of Chicken Entrez Genes (EG) by domain superfamilies (SF).}


\usage{
org.Gg.egSF <- dRDataLoader(RData='org.Gg.egSF')
}

\description{
An R object that contains domain superfamily information for Chicken Entrez Genes. This data is prepared based on SUPERFAMILY database, which provides SCOP domain architecture assignments to all completely sequenced genomes including eukaryotic genomes. The domain architecture for an Entrez gene is the protein product with the longest length of amino acids. Thus, domain superfamily information for Chicken Entrez gene is a list of domain superfamilies (excluding unknown gap) appearing in its domain architecture. 
}

\value{
    an object of class "GS", a list with following components:
    \itemize{
        \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene sets (i.e. SCOP domain superfamilies). The 4 columns are "setID" (i.e. "SCOP domain identifier"), "name" (i.e. "SCOP domain description"), "namespace" (i.e.  "SCOP domain level") and "distance" (i.e.  "SCOP domain classification").
        }
        \item{\code{gs}: a list of gene sets, each storing gene members thereof. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
    }
}

\references{
Morais et al. (2011) SUPERFAMILY 1.75 including a domain-centric gene ontology method. \emph{Nucleic Acids Res}, 39(Database issue):D427-34.
\cr
Andreeva et al. (2008) Data growth and its impact on the SCOP database: new developments. \emph{Nucleic Acids Res}, 36(Database issue):D419-425
}

\keyword{datasets}

\examples{
org.Gg.egSF <- dRDataLoader(RData='org.Gg.egSF')
names(org.Gg.egSF)
}
