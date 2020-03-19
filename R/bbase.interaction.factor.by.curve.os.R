bbase.interaction.factor.by.curve.os <-
function(x, factor, K, bdeg = 3, eps = 1e-5) {
	factor <- droplevels(factor)
	factor.levels <- levels(factor)
	# Ordering in the data set
	ord <- NULL
	for(i in 1:length(factor.levels)) {
		ord <- c(ord, which(factor == factor.levels[i]))
	}
	# Smooth part
	if(length(K) == 1) {
		K <- rep(K, length(factor.levels))
	} else if (length(K) != length(factor.levels)) {
		stop("Error with the number of inner knots for the interaction")
	}
	temp <- interaction.smooth.part <- list()
	for(i in 1:length(factor.levels)) {
		Baux <- bbase.os(x = x[factor == factor.levels[i]], K = K[i], bdeg = bdeg, intercept = FALSE)
		interaction.smooth.part[[i]] <- Baux
		attributes(Baux) <- attributes(Baux)["dim"]
		temp[[i]] <- Baux
	}
	# Join parametric and smooth
	aux <- as.matrix(bdiag(temp))
	B <- aux[order(ord),]
	names(interaction.smooth.part) <- factor.levels
	attr(B,"interaction.smooth.part") <- interaction.smooth.part
	B
}
