\name{Spellman}
\docType{data}
\alias{Spellman}
\title{4 Spellman cancer studies}
\description{
Yeast cell-cycle data set was divided into four subsets which correspond to the four different synchronization methods: alpha arrest (alpha), arrest of cdc15 or cdc28 temperature-sensitive mutant (cdc15 and cdc28), and elutriation (elu). We filtered out genes which have overall missing values >= 10\% or log2 transformed standard deviation >= .45. 1025 genes were left, and the number of time points in the experiments were 18, 24, 17, and 14 for alpha, cdc15, cdc28, and elu, respectively. Additionally, we have imputed missing values using knn. 
}
\usage{Spellman}
\format{A list containing 4 matrices. Each matrix is gene expression data after gene filtering.}
\source{Gene Expression Omnibus (GEO)}
\references{
Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R., Anders, K., Eisen, M. B., Brown, P. O., Botstein, D. and Futcher, B. (1998). Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hybridization. Molecular biology of the cell 9 3273.

}
\keyword{MetaPCA}
