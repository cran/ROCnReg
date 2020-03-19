compute.threshold.FPF.sp <-
function(object, newdata, FPF = 0.5) {
	ncov <- nrow(newdata)
	np <- length(FPF)

	thresholds <- matrix(0, nrow = np, ncol = ncov)
	rownames(thresholds) <- FPF
	fit.new <- predict(object$fit, newdata = newdata)

	if(object$est.cdf == "normal") {
		csf0_inv <- qnorm(1-FPF)
	} else {
		h.residuals <- object$fit$residuals/summary(object$fit)$sigma
		#csf0 <- apply(outer(h.residuals, h.residuals, ">="), 2, mean)
		#csf0_inv <- apply(outer(csf0, FPF, "<="), 2, function(x, z) {
		#	res <- min(c(z[x], max(z)))
		#	res
		#}, z = h.residuals)
		#csf0_inv <- replace(csf0_inv, is.infinite(csf0_inv), max(h.residuals))
		csf0_inv <- quantile(h.residuals, 1-FPF, type = 1)
	}
	for(i in 1:ncov) {
		thresholds[,i] <- fit.new[i] + summary(object$fit)$sigma*csf0_inv
	}
	res <- list()
	res$thresholds <- thresholds
	res
}
