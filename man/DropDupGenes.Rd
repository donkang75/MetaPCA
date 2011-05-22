\name{DropDupGenes}
\alias{DropDupGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MetaPCA: Meta-analysis in the Dimension Reduction of Genomic data
}
\description{
When multiple probesets share the same gene symbols, select only the best probeset in terms of IQR
}
\usage{
DropDupGenes(dat, isParallel=FALSE, nCores=NULL, na.rm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A gene expression matrix which has genes in rows and samples in columns. 
}
  \item{isParallel}{
Whether to use multiple cores in parallel for fast computing. By default, it is false. 
}
  \item{nCores}{
When isParallel is true, the number of cores can be set. By default, all cores in the machine are used in the unix-like machine, and 3 cores are used in windows.
}
  \item{na.rm}{
Whether to remove genes which have no annotation. Default is TRUE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A gene expression matrix which has unique genes in rows and samples in columns. 
}
\references{
Dongwan D. Kang and George C. Tseng. (2011) Meta-PCA: Meta-analysis in the Dimension Reduction of Genomic data. 
}
\author{
Don Kang (donkang75@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	\dontrun{
	#One of example that shows how to generate a expression matrix used in the analysis
	requireAll(c('GEOquery', 'affy', 'hgu133plus2.db'))
	#It might be needed to download the source files first, and save it to local directory
	#such as "./data/Prostate/Varambally" in this example
	#ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE3325/GSE3325_series_matrix.txt.gz
	Varambally <- getGEO('GSE3325', destdir="./data/Prostate/Varambally")
	Varambally <- Varambally[[1]]
	Varambally.sLabel <- as.character(pData(Varambally)$title)
	Varambally.sLabel[grep("Benign",Varambally.sLabel)] <- "Normal"
	Varambally.sLabel[grep("primary",Varambally.sLabel)] <- "Primary"
	Varambally.sLabel[grep("Metastatic",Varambally.sLabel)] <- "Metastasis"
	Varambally <- exprs(Varambally)
	colnames(Varambally) <- Varambally.sLabel
	rownames(Varambally) <- unlist(mget(rownames(Varambally), hgu133plus2SYMBOL))
	Varambally <- DropDupGenes(Varambally, na.rm=TRUE)
	Varambally <- log2(Varambally)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MetaAnalysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
