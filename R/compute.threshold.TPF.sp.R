compute.threshold.TPF.sp <-
function(object, newdata, TPF = 0.5) {
	ncov <- nrow(newdata)
	np <- length(TPF)

	thresholds <- matrix(0, nrow = np, ncol = ncov)
	rownames(thresholds) <- TPF
	fit.new <- predict(object$fit, newdata = newdata)

	if(object$est.cdf == "normal") {
		csf1_inv <- qnorm(1-TPF)
	} else {
		d.residuals <- object$fit$residuals/summary(object$fit)$sigma
		#csf1 <- apply(outer(d.residuals, d.residuals, ">="), 2, mean)
		#csf1_inv <- apply(outer(csf1, TPF, "<="), 2, function(x, z) {
		#	res <- min(c(z[x], max(z)))
		#	res
		#}, z = d.residuals)
		#csf1_inv <- replace(csf1_inv, is.infinite(csf1_inv), max(d.residuals))
		csf1_inv <- quantile(d.residuals, 1-TPF, type = 1)
	}
	for(i in 1:ncov) {
		thresholds[,i] <- fit.new[i] + summary(object$fit)$sigma*csf1_inv
	}
	res <- list()
	res$thresholds <- thresholds
	res
}
