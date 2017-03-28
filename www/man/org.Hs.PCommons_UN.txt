\name{org.Hs.PCommons_UN}
\alias{org.Hs.PCommons_UN}

\title{Human protein interaction undirected network from Pathway Commons (version 7).}

\usage{
org.Hs.PCommons_UN <- dRDataLoader(RData='org.Hs.PCommons_UN')
}

\description{
An igraph object that contains a human protein interaction (undirected) network. The network is extracted from Pathway Commons (version 7; a comprehensive collection of public biological interactions/pathways from multiple data sources). For each interaction and each relation type, data sources in support are: \code{BIND} (see \url{http://bond.unleashedinformatics.com} and \url{http://www.ncbi.nlm.nih.gov/pubmed/12519993}), \code{BioGRID} (see \url{http://thebiogrid.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/16381927}), \code{CORUM} (see \url{http://mips.helmholtz-muenchen.de/genre/proj/corum} and \url{http://www.ncbi.nlm.nih.gov/pubmed/19884131}), \code{DIP} (see \url{http://dip.doe-mbi.ucla.edu} and \url{http://www.ncbi.nlm.nih.gov/pubmed/14681454}), \code{HPRD} (see \url{http://www.hprd.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/18988627}), \code{IntAct} (see \url{http://www.ebi.ac.uk/intact} and \url{http://www.ncbi.nlm.nih.gov/pubmed/24234451}), \code{CTD} (see \url{http://ctdbase.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/25326323}), \code{HumanCyc} (see \url{http://www.humancyc.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/15642094}), \code{KEGG} (see \url{http://www.genome.jp/kegg} and \url{http://www.ncbi.nlm.nih.gov/pubmed/23433509}), \code{PANTHER} (see \url{http://www.pantherdb.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/23193289}), \code{PhosphoSite} (see \url{http://www.phosphosite.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/22135298}), \code{PID} (see \url{http://pid.nci.nih.gov} and \url{http://www.ncbi.nlm.nih.gov/pubmed/18832364}), \code{Reactome} (see \url{http://www.reactome.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/24243840}), \code{ReconX} (see \url{http://humanmetabolism.org} and \url{http://www.ncbi.nlm.nih.gov/pubmed/23455439}), \code{TRANSFAC} (see \url{http://www.gene-regulation.com} and \url{http://www.ncbi.nlm.nih.gov/pubmed/18436575}). Only those interactions with undirected relations (including relation types: "interacts_with", "in_complex_with" and "neighbor_of") are retained. For directed network, please refer to \code{\link{org.Hs.PCommons_DN}}.
}

\value{
    an object of class "igraph" (see \url{http://igraph.org/r/doc/aaa-igraph-package.html}). It has attributes for both vertices and edges. 
    Below are attributes for the vertices:
    \itemize{
        \item{\code{name}: gene symbol for the vertices}
        \item{\code{geneid}: Entrez geneid for the vertices}
        \item{\code{symbol}: gene symbol for the vertices}
        \item{\code{description}: gene description for the vertices}        
    }
    Below are attributes for the edges:
    \itemize{
        \item{\code{in_complex_with}: an integer code indicating support evidence for "in_complex_with" relation, that is, both proteins are members of the same complex}
        \item{\code{interacts_with}: an integer code indicating support evidence for "interacts_with" relation, that is, both proteins are participants of the same interaction}
        \item{\code{neighbor_of}: an integer code indicating support evidence for "neighbor_of" relation, that is, both proteins are participants or controlers of the same interaction}
        \item{Notably, there are two supporting evidences: with the PubMed references (assigned with 100), and the number of different data sources. The integer code is the sum of these two. A value of 0 means no supporting evidences. For examples, a value of 105 means that the relation is supported by the PubMed references and by 5 different data sources.}
    }
}

\references{
Cerami et al. (2011) Pathway Commons, a web resource for biological pathway data. \emph{Nucleic Acids Res}, 39(Database issue):D685-90.
}

\keyword{datasets}

\examples{
org.Hs.PCommons_UN <- dRDataLoader(RData='org.Hs.PCommons_UN')
org.Hs.PCommons_UN
}
