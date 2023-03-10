compute.threshold.TPF.kernel <-
function(object, newdata, TPF = 0.5) {
	ncov <- length(newdata)
	np <- length(TPF)

	thresholds <- matrix(0, nrow = np, ncol = ncov)
	rownames(thresholds) <- TPF
	colnames(thresholds) <- newdata

	fit.mean.new <- npreg(object$bw.mean, exdat = newdata, residuals = TRUE)
	fit.var.new <- npreg(object$bw.var, exdat = newdata, residuals = TRUE)
	d.residuals <- object$fit.mean$resid/sqrt(object$fit.var$mean)

	#csf1 <- apply(outer(d.residuals, d.residuals, ">="), 2, mean)
	#csf1_inv <- apply(outer(csf1, TPF, "<="), 2, function(x, z) {
	#	res <- min(c(z[x], max(z)))
	#	res
	#}, z = d.residuals)
	#csf1_inv <- replace(csf1_inv, is.infinite(csf1_inv), max(d.residuals))

	csf1_inv <- quantile(d.residuals, 1-TPF, type = 1)
	for(i in 1:ncov) {
		thresholds[,i] <- fit.mean.new$mean[i] + sqrt(fit.var.new$mean[i])*csf1_inv
	}
	res <- list()
	res$thresholds <- thresholds
	res
}
