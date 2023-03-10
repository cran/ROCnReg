compute.threshold.TPF.cROC.kernel <-
function(object, newdata, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "cROC.kernel") {
		stop(paste0("This function cannot be used for this object class: ", class(object)[1]))
	}
	# Newdata
	names.cov <- object$covariate

	if(!missing(newdata) && !inherits(newdata, "data.frame"))
		stop("Newdata must be a data frame")
	if(!missing(newdata) && length(names.cov) != 0 &&  sum(is.na(match(names.cov, names(newdata)))))
		stop("Not all needed variables are supplied in newdata") 

	if(missing(newdata)) {
		newdata <- cROCData(object$data, names.cov, object$group)
	} else {
		newdata <- as.data.frame(newdata)
		newdata <- na.omit(newdata[,names.cov,drop = FALSE])
	}

	xp <- newdata[,names.cov]

	res.aux <- compute.threshold.TPF.kernel(object = object$fit$d, newdata = xp, TPF = TPF)
	# Organised results as desired
	thresholds <- vector("list", length(TPF))
	names(thresholds) <- TPF
	for(i in 1:length(TPF)){
		thresholds[[i]] <- matrix(res.aux$thresholds[i,], ncol = 1)
		colnames(thresholds[[i]]) <- "est"
	}
	res <- list()
	res$thresholds <- thresholds
	res$TPF <- TPF
	# Compute associated FPF
	fit.mean.new <- npreg(object$fit$h$bw.mean, exdat = xp, residuals = TRUE)
	fit.var.new <- npreg(object$fit$h$bw.var, exdat = xp, residuals = TRUE)
	h.residuals <- object$fit$h$fit.mean$resid/sqrt(object$fit$h$fit.var$mean)

	aux <- t(t(res.aux$thresholds) - fit.mean.new$mean)
	aux <- t(t(aux)/sqrt(fit.var.new$mean))
	FPF.aux <- matrix(1 - ecdf(h.residuals)(aux), nrow = length(TPF))
	
	# Organised results as desired
	FPF <- vector("list", length(TPF))
	names(FPF) <- TPF
	for(i in 1:length(TPF)){
		FPF[[i]] <- matrix(FPF.aux[i,], ncol = 1)
		colnames(FPF[[i]]) <- "est"
	}
	res$FPF <- FPF
	res$newdata <- newdata
	res 

}
