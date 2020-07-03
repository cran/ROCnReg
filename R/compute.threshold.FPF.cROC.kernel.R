compute.threshold.FPF.cROC.kernel <-
function(object, newdata, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
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
		newdata <- na.omit(newdata[,names.cov,drop = FALSE])
	}

	xp <- newdata[,names.cov]

	res.aux <- compute.threshold.FPF.kernel(object = object$fit$h, newdata = xp, FPF = FPF)
	# Organised results as desired
	thresholds <- vector("list", length(FPF))
	names(thresholds) <- FPF
	for(i in 1:length(FPF)){
		thresholds[[i]] <- matrix(res.aux$thresholds[i,], ncol = 1)
		colnames(thresholds[[i]]) <- "est"
	}
	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	# Compute associated TPF
	fit.mean.new <- npreg(object$fit$d$bw.mean, exdat = xp, residuals = TRUE)
	fit.var.new <- npreg(object$fit$d$bw.var, exdat = xp, residuals = TRUE)
	d.residuals <- object$fit$d$fit.mean$resid/sqrt(object$fit$d$fit.var$mean)

	aux <- t(t(res.aux$thresholds) - fit.mean.new$mean)
	aux <- t(t(aux)/sqrt(fit.var.new$mean))
	TPF.aux <- matrix(1 - ecdf(d.residuals)(aux), nrow = length(FPF))
	
	# Organised results as desired
	TPF <- vector("list", length(FPF))
	names(TPF) <- FPF
	for(i in 1:length(FPF)){
		TPF[[i]] <- matrix(TPF.aux[i,], ncol = 1)
		colnames(TPF[[i]]) <- "est"
	}
	res$TPF <- TPF
	res$newdata <- newdata
	res 

}
