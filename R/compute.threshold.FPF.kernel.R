compute.threshold.FPF.kernel <-
function(object, newdata, FPF = 0.5) {
	ncov <- length(newdata)
	np <- length(FPF)

	thresholds <- matrix(0, nrow = np, ncol = ncov)
	rownames(thresholds) <- FPF
	colnames(thresholds) <- newdata

	fit.mean.new <- npreg(object$bw.mean, exdat = newdata, residuals = TRUE)
	fit.var.new <- npreg(object$bw.var, exdat = newdata, residuals = TRUE)
	h.residuals <- object$fit.mean$resid/sqrt(object$fit.var$mean)

	#csf0 <- apply(outer(h.residuals, h.residuals, ">="), 2, mean)
	#csf0_inv <- apply(outer(csf0, FPF, "<="), 2, function(x, z) {
	#	res <- min(c(z[x], max(z)))
	#	res
	#}, z = h.residuals)
	#csf0_inv <- replace(csf0_inv, is.infinite(csf0_inv), max(h.residuals))

	csf0_inv <- quantile(h.residuals, 1-FPF, type = 1)
	for(i in 1:ncov) {
		thresholds[,i] <- fit.mean.new$mean[i] + sqrt(fit.var.new$mean[i])*csf0_inv
	}
	res <- list()
	res$thresholds <- thresholds
	res
}
