\name{PlotPC2D}
\alias{PlotPC2D}

\title{
MetaPCA: Meta-analysis in the Dimension Reduction of Genomic data
}

\description{
2D PCA plots.
}

\usage{
PlotPC2D(coord, drawObjects=TRUE, drawEllipse=FALSE, dataset.name=NULL, 
		pctInfo=NULL, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, newPlot=TRUE, 
		.points.size=1, .class=rownames(coord), .class.order=NULL, .class.color=NULL, 
		.class2=NULL, .class2.order=NULL, .class2.shape=NULL, .annotation=TRUE, 
		.legend=c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))
}

\arguments{
  \item{coord}{
2D Coordinates matrix of objects. Rows are objects and columns are coordinates.
}
  \item{drawObjects}{
Whether to draw objects as points.
}
  \item{drawEllipse}{
Whether to draw ellipses estimated from objects 2D distribution.
}
  \item{dataset.name}{
Name to be displayed as a part of title.
}
  \item{pctInfo}{
Explained percentage of variance by each PC.
}
  \item{main}{
Main title.
}
  \item{sub}{
Sub title.
}
  \item{xlab}{
Label for x-axis.
}
  \item{ylab}{
Label for y-axis.
}
  \item{newPlot}{
Whether to draw a plot in the new frame.
}
  \item{.points.size}{
Size of objects' points.
}
  \item{.class}{
Object's class label such as disease classification.
}
  \item{.class.order}{
The order of class representation.
}
  \item{.class.color}{
The color of class representation.
}
  \item{.class2}{
The second class label of each object such as study name.
}
  \item{.class2.order}{
The order of 2nd class representation.
}
  \item{.class2.shape}{
The shape of 2nd class representation.
}
  \item{.annotation}{
Whether to present annotation such as x,y axis labels, legend, or titles.
}
  \item{.legend}{
Location of legend in a plot.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NA. A PCA plot is drawn.
}
\references{
Dongwan D. Kang and George C. Tseng. (2011) Meta-PCA: Meta-analysis in the Dimension Reduction of Genomic data. 
}
\author{
Don Kang (donkang75@gmail.com) and George Tseng (ctseng@pitt.edu)
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
	#4 prostate cancer data which have three classes: normal, primary, metastasis
	data(prostate)
	metaPC <- MetaPCA(prostate, method="Angle", doPreprocess=FALSE)
	#Plotting 4 data in the same space with ellipses overlayed!
	coord <- foreach(dd=iter(metaPC$x), .combine=rbind) \%do\% dd$coord
	PlotPC2D(coord[,1:2], drawEllipse=T, dataset.name="Prostate", .class.order=c("Metastasis","Primary","Normal"), 
			.class.color=c('red','#838383','blue'), .annotation=T, newPlot=T,
			.class2=rep(names(metaPC$x), times=sapply(metaPC$x,function(x)nrow(x$coord))), 
			.class2.order=names(metaPC$x), .points.size=1)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PCA }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
