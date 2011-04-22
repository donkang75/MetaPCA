PreprocessMetaAnalysis <- function(DList, cutRatioByMean=.4, cutRatioByVar=.4, doImpute=FALSE, na.rm.pct=.1, na.rm.pct.each=.5, verbose=FALSE) {	
	studies <- names(DList)
	
	DList <- foreach(dat=iter(DList)) %do% {
		dat[rowSums(is.na(dat))<=ncol(dat)*na.rm.pct.each,] #remove if more than some pct is missing in single study
	}
	
	.genes <- foreach(dat=iter(DList)) %do% {rownames(dat)}
	intersect.rec <- function(.list, ...){
		if(length(.list)==1) return(.list[[1]])
		Recall(c(list(intersect(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
	}
	.genes <- na.omit(intersect.rec(.genes))
	
	printLog(paste("# of matched genes :", length(.genes)), verbose)
	
	DD <- foreach(dat=iter(DList), .combine=cbind) %do% {
		dat[match(.genes,rownames(dat)),]
	}
	
	DD <- DD[rowSums(is.na(DD))<ncol(DD)*na.rm.pct,]
	
	n <- c(0, cumsum(sapply(DList, ncol)))
	for(i in 1:(length(n)-1)) {
		DList[[i]] <- DD[,(n[i]+1):n[i+1]]
	}
	
	printLog(paste("# of genes after na.rm threshold:", nrow(DList[[1]])), verbose)
	
	if(doImpute) {
		requireAll("impute")
		for(i in 1:length(studies)) {
			if(any(is.na(DList[[i]])))
				DList[[i]] <- impute.knn(DList[[i]])$data
		}
	}
	
	numLeft <- floor(nrow(DList[[1]])*(1-cutRatioByMean)) 
	
	#rank sum of row means
	DList.rankM <- as.matrix(foreach(dat=iter(DList), .combine=cbind) %do% {
				rank(rowMeans(dat))
			})
	DList <- foreach(dat=iter(DList)) %do% {
		dat[order(rowSums(DList.rankM),decreasing=TRUE)[1:numLeft],]
	} 
	
	numLeft <- floor(numLeft*(1-cutRatioByVar))
	
	#rank sum of row var
	DList.rankV <- as.matrix(foreach(dat=iter(DList), .combine=cbind) %do% {
				rank(rowVars(dat))
			})
	DList <- foreach(dat=iter(DList)) %do% {
		dat[order(rowSums(DList.rankV),decreasing=TRUE)[1:numLeft],]
	} 
	
	names(DList) <- studies
	return(DList)
}

#class: usually disease labels
#class2: different sources like multiple studies
PlotPC2D <- function(coord, drawObjects=TRUE, drawEllipse=FALSE, dataset.name=NULL, 
		pctInfo=NULL, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, newPlot=TRUE, .points.size=1,
		.class=rownames(coord), .class.order=NULL, .class.color=NULL, .class2=NULL, .class2.order=NULL, .class2.shape=NULL, 
		.annotation=TRUE, .legend=c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")) 
{
	.legend <- match.arg(.legend)
	if(is.null(.class.order)) {
		.label <- factor(.class) 		
	} else {
		.label <- factor(.class, levels=.class.order)
	}
	if(is.null(.class.color)) {
		.class.color <- 1:nlevels(.label) 		
	} 
	
	if(!is.null(.class2)) {
		if(is.null(.class2.order)) {
			.label2 <- factor(.class2) 		
		} else {
			.label2 <- factor(.class2, levels=.class2.order)
		}
		if(is.null(.class2.shape)) {
			.class2.shape <- 1:nlevels(.label2)  		
		} 
	}
	
	colnames(coord) <- c("PC1","PC2")
	
	if(newPlot) {
		.sysname <- Sys.info()["sysname"]
		if(.sysname=="Windows")
			windows()
		else if(.sysname=="Darwin")
			quartz()
		else #*nix
			X11()
	}
	
	
	if(drawEllipse) {
		requireAll(c("ellipse","MASS"))
		.ellip <- NULL 
		for(lev in levels(.label)) {
			.mve <- cov.mve(coord[.label==lev,1:2], method="classical")
			.ellip <- rbind(.ellip, ellipse(x=.mve$cov, centre=.mve$center, npoints=100))
		}
		plot(rbind(.ellip,coord), type="n", xlab="", ylab="", axes=FALSE, xlim=range(rbind(.ellip,coord)), ylim=range(rbind(.ellip,coord)))		
		for(i in 1:nlevels(.label)) {
			lines(.ellip[(1+100*(i-1)):(100*i),], col=.class.color[i])
		}
	} else {
		plot(coord, type="n", xlab="", ylab="", axes=FALSE, xlim=range(coord), ylim=range(coord))	
	}
	
	if(drawObjects) 
		for(i in 1:nlevels(.label)) {
			if(exists('.label2')) {
				for(j in 1:nlevels(.label2)) 
					points(coord[.label==levels(.label)[i] & .label2==levels(.label2)[j],], pch=.class2.shape[j], col=.class.color[i], cex=.points.size)
			} else
				points(coord[.label==levels(.label)[i],], pch=20, col=.class.color[i], cex=.points.size)			
		}
	
	if(.annotation) {
		axis(side=1,lwd=4,tck=-0.02)
		axis(side=2,lwd=4,tck=-0.02)
	}
	box(bty="L",lwd=4)
	abline(v=0, lty=2)
	abline(h=0, lty=2)	
	
	if(.annotation) {
		if(exists('.label2')) {
			legend2(.legend, legend=c(levels(.label),levels(.label2)), pch=c(rep(-1,nlevels(.label)), .class2.shape), fill=c(.class.color, rep(-1,nlevels(.label2))))
		} else
			legend(.legend, legend=levels(.label), text.col=.class.color, pch=20, col=.class.color)
		
		.pctInfo <- c(ifelse(is.null(pctInfo), "", paste(" (",pctInfo[1],"%)",sep="")),
				ifelse(is.null(pctInfo), "", paste(" (",pctInfo[2],"%)",sep="")))
		
		.main <- ifelse(is.null(main), paste("Meta-PCA [", dataset.name, "]", sep=""), main)
		.xlab <- ifelse(is.null(xlab), paste("Meta-PC1", .pctInfo[1], sep=""), xlab)
		.ylab <- ifelse(is.null(ylab), paste("Meta-PC2", .pctInfo[2], sep=""), ylab)
		
		title(main=.main, sub=sub, xlab=.xlab, ylab=.ylab)
	}
	
}

#from matrixStats
colSds <- function (x, center = NULL, ...) {
	x <- rowVars(t(x), ...)
	sqrt(x)
}

rowVars <- function (x, center = NULL, ...) 
{
	n <- !is.na(x)
	n <- rowSums(n)
	n[n <= 1] <- NA
	if (is.null(center)) {
		center <- rowMeans(x, ...)
	}
	x <- x - center
	x <- x * x
	x <- rowSums(x, ...)
	x <- x/(n - 1)
	x
}

#maxFeatures at most
MetaSoftThreshold <- function(aMat, maxFeatures, lambda){
	.rs <- rowSums(abs(aMat))
	
	if(!is.null(maxFeatures)) 
		lambda <- max(lambda, sort(.rs, decreasing=TRUE)[maxFeatures+1])
	
	aMat[which(.rs <= lambda),] <- 0 #less than overall threshold
	
	tmp <- abs(aMat) - lambda/ncol(aMat) #substract individual threshold
	tmp[tmp<0] <- 0 #less than individual threshold
	sign(aMat) * tmp 
}

#pxn 
NormalizeMatrix <- function(mat) {
	.norm <- sqrt(colSums(mat^2))
	.norm <- ifelse(.norm==0, 1, .norm)
	sweep(mat, 2, .norm, '/')
}


legend2 <- function (x, y = NULL, legend, fill = NULL, col = par("col"), 
		border = "black", lty, lwd, pch, angle = 45, density = NULL, 
		bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
		box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, 
		xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 
				0.5), text.width = NULL, text.col = par("col"), merge = do.lines && 
				has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE, 
		title = NULL, inset = 0, xpd, title.col = text.col, title.adj = 0.5, 
		seg.len = 2) 
{
	if (missing(legend) && !missing(y) && (is.character(y) || 
				is.expression(y))) {
		legend <- y
		y <- NULL
	}
	mfill <- !missing(fill) || !missing(density)
	if (!missing(xpd)) {
		op <- par("xpd")
		on.exit(par(xpd = op))
		par(xpd = xpd)
	}
	title <- as.graphicsAnnot(title)
	if (length(title) > 1) 
		stop("invalid title")
	legend <- as.graphicsAnnot(legend)
	n.leg <- if (is.call(legend)) 
				1
			else length(legend)
	if (n.leg == 0) 
		stop("'legend' is of length 0")
	auto <- if (is.character(x)) 
				match.arg(x, c("bottomright", "bottom", "bottomleft", 
								"left", "topleft", "top", "topright", "right", "center"))
			else NA
	if (is.na(auto)) {
		xy <- xy.coords(x, y)
		x <- xy$x
		y <- xy$y
		nx <- length(x)
		if (nx < 1 || nx > 2) 
			stop("invalid coordinate lengths")
	}
	else nx <- 0
	xlog <- par("xlog")
	ylog <- par("ylog")
	rect2 <- function(left, top, dx, dy, density = NULL, angle, 
			...) {
		r <- left + dx
		if (xlog) {
			left <- 10^left
			r <- 10^r
		}
		b <- top - dy
		if (ylog) {
			top <- 10^top
			b <- 10^b
		}
		rect(left, top, r, b, angle = angle, density = density, 
				...)
	}
	segments2 <- function(x1, y1, dx, dy, ...) {
		x2 <- x1 + dx
		if (xlog) {
			x1 <- 10^x1
			x2 <- 10^x2
		}
		y2 <- y1 + dy
		if (ylog) {
			y1 <- 10^y1
			y2 <- 10^y2
		}
		segments(x1, y1, x2, y2, ...)
	}
	points2 <- function(x, y, ...) {
		if (xlog) 
			x <- 10^x
		if (ylog) 
			y <- 10^y
		points(x, y, ...)
	}
	text2 <- function(x, y, ...) {
		if (xlog) 
			x <- 10^x
		if (ylog) 
			y <- 10^y
		text(x, y, ...)
	}
	if (trace) 
		catn <- function(...) do.call("cat", c(lapply(list(...), 
									formatC), list("\n")))
	cin <- par("cin")
	Cex <- cex * par("cex")
	if (is.null(text.width)) 
		text.width <- max(abs(strwidth(legend, units = "user", 
								cex = cex)))
	else if (!is.numeric(text.width) || text.width < 0) 
		stop("'text.width' must be numeric, >= 0")
	xc <- Cex * xinch(cin[1L], warn.log = FALSE)
	yc <- Cex * yinch(cin[2L], warn.log = FALSE)
	if (xc < 0) 
		text.width <- -text.width
	xchar <- xc
	xextra <- 0
	yextra <- yc * (y.intersp - 1)
	ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
	ychar <- yextra + ymax
	if (trace) 
		catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
						ychar))
	if (mfill) {
		xbox <- xc * 0.8
		ybox <- yc * 0.5
		dx.fill <- xbox
	}
	do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
									0))) || !missing(lwd)
	n.legpercol <- if (horiz) {
				if (ncol != 1) 
					warning("horizontal specification overrides: Number of columns := ", 
							n.leg)
				ncol <- n.leg
				1
			}
			else ceiling(n.leg/ncol)
	has.pch <- !missing(pch) && length(pch) > 0
	if (do.lines) {
		x.off <- if (merge) 
					-0.7
				else 0
	}
	else if (merge) 
		warning("'merge = TRUE' has no effect when no line segments are drawn")
	if (has.pch) {
		if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
				type = "c") > 1) {
			if (length(pch) > 1) 
				warning("not using pch[2..] since pch[1L] has multiple chars")
			np <- nchar(pch[1L], type = "c")
			pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
		}
	}
	if (is.na(auto)) {
		if (xlog) 
			x <- log10(x)
		if (ylog) 
			y <- log10(y)
	}
	if (nx == 2) {
		x <- sort(x)
		y <- sort(y)
		left <- x[1L]
		top <- y[2L]
		w <- diff(x)
		h <- diff(y)
		w0 <- w/ncol
		x <- mean(x)
		y <- mean(y)
		if (missing(xjust)) 
			xjust <- 0.5
		if (missing(yjust)) 
			yjust <- 0.5
	}
	else {
		h <- (n.legpercol + (!is.null(title))) * ychar + yc
		w0 <- text.width + (x.intersp + 1) * xchar
		if (mfill) 
			w0 <- w0 + dx.fill
		if (do.lines) 
			w0 <- w0 + (seg.len + +x.off) * xchar
		w <- ncol * w0 + 0.5 * xchar
		if (!is.null(title) && (abs(tw <- strwidth(title, units = "user", 
									cex = cex) + 0.5 * xchar)) > abs(w)) {
			xextra <- (tw - w)/2
			w <- tw
		}
		if (is.na(auto)) {
			left <- x - xjust * w
			top <- y + (1 - yjust) * h
		}
		else {
			usr <- par("usr")
			inset <- rep(inset, length.out = 2)
			insetx <- inset[1L] * (usr[2L] - usr[1L])
			left <- switch(auto, bottomright = , topright = , 
					right = usr[2L] - w - insetx, bottomleft = , 
					left = , topleft = usr[1L] + insetx, bottom = , 
					top = , center = (usr[1L] + usr[2L] - w)/2)
			insety <- inset[2L] * (usr[4L] - usr[3L])
			top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
							h + insety, topleft = , top = , topright = usr[4L] - 
							insety, left = , right = , center = (usr[3L] + 
								usr[4L] + h)/2)
		}
	}
	if (plot && bty != "n") {
		if (trace) 
			catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
					h, ", ...)", sep = "")
		rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
				lwd = box.lwd, lty = box.lty, border = box.col)
	}
	xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1), 
				rep.int(n.legpercol, ncol)))[1L:n.leg]
	yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol, 
						ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
	if (mfill) {
		if (plot) {
			fill <- rep(fill, length.out = n.leg)
			rect2(left = (xt + dx.fill/2)[which(fill!=-1)], top = (yt + ybox/2)[which(fill!=-1)], dx = xbox, dy = ybox, 
					col = fill[which(fill!=-1)], density = density, angle = angle, 
					border = border)
		}
		xt <- xt + dx.fill
	}
	if (plot && (has.pch || do.lines)) 
		col <- rep(col, length.out = n.leg)
	if (missing(lwd)) 
		lwd <- par("lwd")
	if (do.lines) {
		if (missing(lty)) 
			lty <- 1
		lty <- rep(lty, length.out = n.leg)
		lwd <- rep(lwd, length.out = n.leg)
		ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
		if (trace) 
			catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
					yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
		if (plot) 
			segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
							xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
					col = col[ok.l])
		xt <- xt + (seg.len + x.off) * xchar
	}
	if (has.pch) {
		pch <- rep(pch, length.out = n.leg)
		pt.bg <- rep(pt.bg, length.out = n.leg)
		pt.cex <- rep(pt.cex, length.out = n.leg)
		pt.lwd <- rep(pt.lwd, length.out = n.leg)
		ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
		x1 <- (if (merge && do.lines) 
				xt - (seg.len/2) * xchar
			else xt)[ok]
		y1 <- yt[ok]
		if (trace) 
			catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
					", ...)")
		if (plot) 
			points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
					bg = pt.bg[ok], lwd = pt.lwd[ok])
	}
	xt <- xt + x.intersp * xchar
	if (plot) {
		if (!is.null(title)) 
			text2(left + w * title.adj, top - ymax, labels = title, 
					adj = c(title.adj, 0), cex = cex, col = title.col)
		text2(xt, yt, labels = legend, adj = adj, cex = cex, 
				col = text.col)
	}
	invisible(list(rect = list(w = w, h = h, left = left, top = top), 
					text = list(x = xt, y = yt)))
}

printLog <- function(msg, verbose) {
	if(verbose)
		print(sprintf("[%s] %s", Sys.time(), msg))
}

l1median_HoCr2 <- function(X, ...) {
	pcaPP:::l1median_HoCr(X, ...)$par
}

l1median_VaZh2 <- function(X, ...) {
	pcaPP:::l1median_VaZh(X, ...)$par
}
