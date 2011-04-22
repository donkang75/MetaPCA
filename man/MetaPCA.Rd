\name{MetaPCA}
\alias{MetaPCA}

\title{
MetaPCA: Meta-analysis in the Dimension Reduction of Microarray data
}

\description{
MetaPCA implements simultaneous dimension reduction using PCA when multiple studies are combined. We propose two basic ideas to find a common PC subspace by eigenvalue maximization approach and angle minimization approach, and we extend the concept to incorporate Robust PCA and Sparse PCA in the meta-analysis realm.
}

\usage{
MetaPCA(DList, method=c("Angle","Eigen","RobustAngle","SparseAngle"), robust.var=c("qn","mad"), nPC=2,
			.weight=rep(1/length(DList),length(DList)), sparse.maxFeatures=NULL, sparse.lambda=NULL, 
			sparse.max.iter=100, sparse.eps=1e-3, .scale=FALSE, .scaleAdjust=TRUE, doPreprocess=TRUE, 
			cutRatioByMean=.4, cutRatioByVar=.4, doImpute=TRUE,	na.rm.pct=.1, na.rm.pct.each=.5, 
			verbose=FALSE)
}

\arguments{
  \item{DList}{
A list of all data matrices; Each data name should be set as the name of each list element. Each data should be a numeric matrix that has genes in the rows and samples in the columns. Row names should be official gene symbols and column names be sample labels.
}
  \item{method}{
A vector of four meta PCA methods. The first two methods are basic approaches; the last two are extended approaches of robust PCA and sparse PCA but may be rather slower than the basic methods. Default is "Angle", which is angle minimization method. See the details in the reference.
}
  \item{robust.var}{
Robust measure of variance when "RobustAngle" method was selected in the method. 
}
  \item{nPC}{
The number of returned PC's, i.e. the number of dimension reduced by PCA.
}
  \item{.weight}{
Weight for each data if information is available. Default is equal weight.
}
  \item{sparse.maxFeatures}{
The number of genes left for the Sparse PCA approach. If NULL (default), it is determined based on the default lambda.
}
  \item{sparse.lambda}{
The parameter lambda which determines the sparsity of loading vectors. The default is calculated as the number of data divided by square root of the number of overall genes.
}
  \item{sparse.max.iter}{
The number of maximum iteration for achieving convergence of sparse loading vectors. Default is 100.
}
  \item{sparse.eps}{
The convergence decision precision level. Default is 1e-3.
}
  \item{.scale}{
Whether to apply gene based normalization. Default is FALSE. But for the "Eigen" method, gene scaling is recommended for the comparability reason of covariance matrix.
}
  \item{.scaleAdjust}{
Whether to apply scaling adjustment for a comparable visualization. Default is TRUE.
}
  \item{doPreprocess}{
Whether to apply gene filtering. Default is TRUE. However "SparseAngle" method do not use gene filtering.
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
list object having the specified number of PC's of all data sets and loading matrix of meta subspace.
}
\references{
Dongwan D. Kang and George C. Tseng. (2011) Meta-PCA: Meta-analysis in the Dimension Reduction of Microarray data. 
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
	#Spellman, 1998 Yeast cell cycle data set
	#Consider each synchronization method as a separate data
	data(Spellman) 
	pc <- list(alpha=prcomp(t(Spellman$alpha))$x, cdc15=prcomp(t(Spellman$cdc15))$x,
			cdc28=prcomp(t(Spellman$cdc28))$x, elu=prcomp(t(Spellman$elu))$x)
	#There are currently 4 meta-pca methods. Run either one of following four.
	metaPC <- MetaPCA(Spellman, method="Eigen", doPreprocess=FALSE)
	metaPC <- MetaPCA(Spellman, method="Angle", doPreprocess=FALSE)
	metaPC <- MetaPCA(Spellman, method="RobustAngle", doPreprocess=FALSE)
	metaPC <- MetaPCA(Spellman, method="SparseAngle", doPreprocess=FALSE)
	#Comparing between usual pca and meta-pca
	#The first lows are four data sets based on usual PCA, and 
	#the second rows are by MetaPCA
	#We're looking for a cyclic pattern.
	par(mfrow=c(2,4), cex=1, mar=c(0.2,0.2,0.2,0.2))
	for(i in 1:4) {
		plot(pc[[i]][,1], pc[[i]][,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n")
		text(pc[[i]][,1], pc[[i]][,2], 1:nrow(pc[[i]]), cex=1.5)
		lines(pc[[i]][,1], pc[[i]][,2])
	}
	for(i in 1:4) {
		plot(metaPC$x[[i]]$coord[,1], metaPC$x[[i]]$coord[,2], type="n", xlab="", ylab="", xaxt="n", yaxt="n")
		text(metaPC$x[[i]]$coord[,1], metaPC$x[[i]]$coord[,2], 1:nrow(metaPC$x[[i]]$coord), cex=1.5)
		lines(metaPC$x[[i]]$coord[,1], metaPC$x[[i]]$coord[,2])
	}

	#4 prostate cancer data which have three classes: normal, primary, metastasis
	data(prostate)
	#There are currently 4 meta-pca methods. Run either one of following four.
	metaPC <- MetaPCA(prostate, method="Eigen", doPreprocess=FALSE, .scale=TRUE)
	metaPC <- MetaPCA(prostate, method="Angle", doPreprocess=FALSE)
	metaPC <- MetaPCA(prostate, method="RobustAngle", doPreprocess=FALSE)
	metaPC <- MetaPCA(prostate, method="SparseAngle", doPreprocess=FALSE)
	#Plotting 4 data in the same space!
	coord <- foreach(dd=iter(metaPC$x), .combine=rbind) %do% dd$coord
	PlotPC2D(coord[,1:2], drawEllipse=F, dataset.name="Prostate", .class.order=c("Metastasis","Primary","Normal"), 
			.class.color=c('red','#838383','blue'), .annotation=T, newPlot=T,
			.class2=rep(names(metaPC$x), times=sapply(metaPC$x,function(x)nrow(x$coord))), 
			.class2.order=names(metaPC$x), .points.size=1)

	#In the case of "SparseAngle" method, the top contributing genes for all studies can be determined
	#For instance, top 20 genes in 1st PC and their coefficients
	metaPC$v[order(abs(metaPC$v[,1]), decreasing=TRUE),1][1:20] 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PCA }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
