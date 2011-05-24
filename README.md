MetaPCA: Meta-analysis in the Dimension Reduction of Genomic data
====================================================================

Introduction
------------
__MetaPCA__ implements simultaneous dimension reduction using Principal Component Analysis (PCA) when multiple studies are combined. We propose two basic ideas to find a common PC subspace by eigenvalue maximization approach and angle minimization approach, and we extend the concept to incorporate Robust PCA and Sparse PCA in the meta-analysis realm.

__MetaPCA__ implements our proposed four methods to find common subspace among multiple studies (datasets): 

* (1) Eigenvalue Maximization Approach
* (2) Angle Minimization Approach
* (3) Robust Angle Minimization Approach which extends Robust PCA
* (4) Sparse Angle Minimization Approach which extends Sparse PCA 

For the detailed information, please see the references.

Installation
--------------
To install this package, save a proper package file for the target OS to the working directory, then run:

### Windows            
[MetaPCA_0.1.3.zip] (http://cran.r-project.org/bin/windows/contrib/r-release/MetaPCA_0.1.3.zip)

        install.packages("MetaPCA_0.1.3.zip", repos=NULL, type="win.binary")

### Mac OS X            
[MetaPCA_0.1.3.tgz] (http://cran.r-project.org/bin/macosx/leopard/contrib/r-release/MetaPCA_0.1.3.tgz)

        install.packages("MetaPCA_0.1.3.tgz", repos=NULL, type="mac.binary")

### Linux            
[MetaPCA_0.1.3.tar.gz] (http://cran.r-project.org/src/contrib/MetaPCA_0.1.3.tar.gz)

        install.packages("MetaPCA_0.1.3.tar.gz", repos=NULL, type="source")

Examples
-------------
        library(MetaPCA)
        requireAll("foreach")

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
		
References
----------
Dongwan D. Kang and George C. Tseng. (2011) Meta-PCA: Meta-analysis in the Dimension Reduction of Genomic data. (In preparation) 

