\name{MetaPCA-package}
\alias{MetaPCA-package}
\docType{package}
\title{
MetaPCA: Meta-analysis in the Dimension Reduction of Microarray data
}
\description{
MetaPCA implements simultaneous dimension reduction using PCA when multiple studies are combined. We propose two basic ideas to find a common PC subspace by eigenvalue maximization approach and angle minimization approach, and we extend the concept to incorporate Robust PCA and Sparse PCA in the meta-analysis realm.
}
\details{
\tabular{ll}{
Package: \tab MetaPCA\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2011-04-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Don Kang <donkang75@gmail.com>
}
\references{
Dongwan D. Kang and George C. Tseng. (2011) Meta-PCA: Meta-analysis in the Dimension Reduction of Microarray data. 
}
\keyword{ PCA }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
\seealso{
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
