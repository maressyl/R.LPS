## Linear Predictor Score plot
## Author : Sylvain Mareschal <maressyl@gmail.com>
plot.LPS <- function(x, y, method=c("Wright", "Radmacher", "exact"), xlim, yaxt="n", xlab="LPS", ylab, ...) {
	# X
	if(missing(xlim)) xlim <- c(min(x$means - 4*x$sds), max(x$means + 4*x$sds))
	xval <- seq(from=xlim[1], to=xlim[2], length.out=1000)
	
	# Y
	method <- match.arg(method)
	if(missing(y)) y <- "density"
	y <- match.arg(y, c("density", "probability"))
	if(missing(ylab)) ylab <- sprintf("%s (%s)", y, method)
	if(y == "density") {
		# Density plot
		if(method == "Wright") {
			# Wright et al (gaussian densities)
			yval1 <- dnorm(xval, mean=x$means[1], sd=x$sds[1])
			yval2 <- dnorm(xval, mean=x$means[2], sd=x$sds[2])
			type <- "l"
		} else {
			# Other (discrete values)
			xval <- c(x$scores[[1]], x$scores[[2]])
			yval1 <- c(rank(x$scores[[1]]), rep(NA, length(x$scores[[2]])))
			yval2 <- c(rep(NA, length(x$scores[[1]])), rank(x$scores[[2]]))
			type <- "p"
		}
	} else {
		# Probability plot
		p <- predict(x, newdata=xval, type="probability", method=method, plot=FALSE)
		yval1 <- p[,1]
		yval2 <- p[,2]
		type <- "l"
	}
	ylim <- range(c(yval1, yval2), na.rm=TRUE)
	
	# Plot
	plot(x=xval, y=yval1, type=type, col="blue", xlim=xlim, ylim=ylim, yaxt=yaxt, xlab=xlab, ylab=ylab, ...)
	par(new=TRUE)
	plot(x=xval, y=yval2, type=type, col="red", xlim=xlim, ylim=ylim, yaxt=yaxt, xlab=xlab, ylab=ylab, ...)
	
	# Radmacher's means
	if(method == "Radmacher") {
		abline(v=x$means[1], col="blue")
		abline(v=x$means[2], col="red")
	}
	
	# Legend
	legend(
		x = "topleft",
		inset = 0.01,
		bg = "#DDDDDD",
		legend = c(x$classes),
		col = c("blue", "red"),
		lty = "solid"
	)
}
