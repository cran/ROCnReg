predictive.checks.helper <-
function(y, yrep, marker, statistics, devnew, group = c("D", "H"), ndensity = NULL) {
	group <- match.arg(group)
	main.stat <- paste0(statistics, " (Group ", group, ")")
	main.hist <- paste0("Density (Group ", group, ")")
	xlab.hist <- paste0(marker, " (Group ", group, ")")

	nrep <- ncol(yrep)

	i <- 1
	for(stat in statistics) {
		if(i != 1 & devnew) dev.new()
		yrepstat <- apply(yrep, 2, function(y, stat) {do.call(stat, list(y))}, stat = stat)
		yrepstat <- trim.out(yrepstat, trim = 0.025)
		ystat <- do.call(stat, list(y))
		xlim <- range(c(yrepstat,ystat))
		hist(yrepstat, col = "gray60", main = main.stat[i], xlim = xlim, xlab = "Statistic")
		abline(v = ystat,col="red",lwd=3)
		i = i + 1
	}
	if(devnew) dev.new()

	if(is.null(ndensity)) ndensity <- 512

	ylim <- c(0, max(density(y)$y) + 0.2*max(density(y)$y))
	xlim <- c(min(density(y)$x) - 0.2, max(density(y)$x) - 0.2)
	plot(density(yrep[,1], n = ndensity), col = "lightskyblue1", ylim = ylim, xlim = xlim, main = main.hist, xlab = xlab.hist)
	# Only a sample
	s <- sample(1:nrep, ifelse(nrep < 500, nrep, 500)) 
	for(i in s){
		lines(density(yrep[,i], n = ndensity), col = "lightskyblue1")  
	}
	lines(density(y),col = "black",lwd = 4)
}
