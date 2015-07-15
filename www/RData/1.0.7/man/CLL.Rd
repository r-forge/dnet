\name{CLL}
\alias{CLL}

\title{Transcriptome dataset for Chronic Lymphocytic Leukemia patients from Chuang et al. (2012)}

\usage{
CLL <- dRDataLoader(RData='CLL')
}

\description{
This dataset involves 130 patients with chronic lymphocytic leukemia (CLL). When enrolled in the study, these CLL patients had not received prior therapy for CLL. Additional covariate about the time to treatment (i.e. prognosis) is available. The dataset has been normalised and log2-transformed, and provided as an 'ExpressionSet' object.
}

\value{
    an object of class "ExpressionSet". It has slots for "assayData", "phenoData", and "featureData": 
    \itemize{
        \item{\code{assayData}: a matrix of 54675 features X 130 samples}
        \item{\code{phenoData}: variables describing sample phenotypes (i.e. columns in assayData), including information about samples: "Name" for sample names, "Time" for sampling time to first treatment (years) and "Treatment" for treatment event (1:yes, 0:no)}
        \item{\code{featureData}: variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol and "Desc" for gene description}
    }
}

\references{
Chuang et al. (2012). Subnetwork-based analysis of chronic lymphocytic leukemia identifies pathways that associate with disease progression. \emph{Blood}, 120(13):2639-49.
}

\keyword{datasets}

\examples{
CLL <- dRDataLoader(RData='CLL')
CLL
library(Biobase)
# extract information about the first 5 samples
pData(CLL)[1:5,]
# extract information about the first 5 features
fData(CLL)[1:5,]
}
