\name{PreprocessMetaAnalysis}
\alias{PreprocessMetaAnalysis}

\title{
MetaPCA: Meta-analysis in the Dimension Reduction of Genomic data
}

\description{
Preprocessing for microarray meta-analysis. It is about gene filtering and missing value imputation.
}

\usage{
PreprocessMetaAnalysis(DList, cutRatioByMean=.4, cutRatioByVar=.4, doImpute=FALSE, na.rm.pct=.1, na.rm.pct.each=.5, verbose=FALSE)
}

\arguments{
  \item{DList}{
A list of all data matrices; Each data name should be set as the name of each list element. Each data should be a numeric matrix that has genes in the rows and samples in the columns. Row names should be official gene symbols and column names be sample labels.
}
  \item{cutRatioByMean}{
Proportion of genes filtered by study-wise mean. Default is 40\%.
}
  \item{cutRatioByVar}{
Proportion of genes filtered by study-wise variance. Default is 40\%.
}
  \item{doImpute}{
Whether to impute missing genes. Default is TRUE, and default imputation method is knn.
}
  \item{na.rm.pct}{
Proportion of genes filtered by study-wise missing proportion. Default is 10\%.
}
  \item{na.rm.pct.each}{
Proportion of genes filtered by each study's missing proportion. Default is 50\%.
}
  \item{verbose}{
Whether to print logs. Default is FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list object of all data matrices after filtering and imputation.
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
DList <- PreprocessMetaAnalysis(list(Yu=Yu, Lapointe=Lapointe, Tomlins=Tomlins, Varambally=Varambally), 
		cutRatioByMean=.1, cutRatioByVar=.1, doImpute=T, na.rm.pct=.2) #3056
str(DList)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PCA }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
