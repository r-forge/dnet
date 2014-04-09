\name{TCGA_mutations}
\alias{TCGA_mutations}

\title{TCGA mutational profiles across 12 major cancer types from Kandoth et al. (2013)}

\usage{
load(url("http://dnet.r-forge.r-project.org/data/Datasets/TCGA_mutations.RData"))
}

\description{
This dataset is available from TCGA, containing somatic mutational profiles for 3130 cancer samples (i.e. only those having clinical data). These cancer samples belong to one of 12 major cancer types, including breast adenocarcinoma (BRCA), lung adenocarcinoma (LUAD), lung squamous cell carcinoma (LUSC), uterine corpus endometrial carcinoma (UCEC), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC), colon and rectal carcinoma (COAD/READ), bladder urothelial carcinoma (BLCA), kidney renal clear cell carcinoma (KIRC), ovarian serous carcinoma (OV) and acute myeloid leukaemia (LAML). For each patient sample, somatic mutations are represented as a profile of binary (1, 0) states on genes, where '1' indicates a gene for which mutation has occurred in the tumor relative to germ line. The dataset is provided as an 'ExpressionSet' object.
}

\value{
    an object of class "ExpressionSet". It has slots for "assayData", "phenoData", and "featureData": 
    \itemize{
        \item{\code{assayData}: a matrix of 19477 genes X 3130 samples}
        \item{\code{phenoData}: variables describing sample phenotypes (i.e. columns in assayData), including information about samples: "Age", "Days_to_death", "Days_to_last_followup", "Vita_status", "Gender", "Date_diagnosis", "TCGA_tumor_type", "Tumor_stage", "Tumor_grade"}
        \item{\code{featureData}: variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol, "Desc" for gene description, "Synonyms" for gene symbol alias}
    }
}

\references{
Kandoth et al. (2013). Mutational landscape and significance across 12 major cancer types. \emph{Nature}, 502(7471):333-9.
}

\keyword{datasets}
\examples{
load(url("http://dnet.r-forge.r-project.org/data/Datasets/TCGA_mutations.RData"))
TCGA_mutations
# extract information about the first 5 samples
pData(TCGA_mutations)[1:5,]
# extract information about the first 5 features
fData(TCGA_mutations)[1:5,]
# number of samples for each cancer type
table(pData(TCGA_mutations)$TCGA_tumor_type)
}
