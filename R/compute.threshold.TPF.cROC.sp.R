compute.threshold.TPF.cROC.sp <-
function(object, newdata, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "cROC.sp") {
		stop(paste0("This function cannot be used for this object class: ", class(object)[1]))
	}

	# Newdata
	names.cov.h <- all.vars(object$formula$h)[-1]
	names.cov.d <- all.vars(object$formula$d)[-1]
	names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])

	if(!missing(newdata) && !inherits(newdata, "data.frame"))
		stop("Newdata must be a data frame")
	if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
		stop("Not all needed variables are supplied in newdata") 

	if(missing(newdata)) {
		newdata <- cROCData(object$data, names.cov, object$group)
	} else {
		newdata <- as.data.frame(newdata)
		newdata <- na.omit(newdata[,names.cov,drop = FALSE])
	}  
	res.aux <- compute.threshold.TPF.sp(object = list(est.cdf = object$est.cdf, fit = object$fit$d), newdata = newdata, TPF = TPF)
	
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
	
	fit.new <- predict(object$fit$h, newdata = newdata)
	aux <- t(t(res.aux$thresholds) - fit.new)
	aux <- t(t(aux)/summary(object$fit$h)$sigma)

	if(object$est.cdf == "normal") {
		FPF.aux <- 1 - pnorm(aux)
	} else {
		h.residuals <- object$fit$h$residuals/summary(object$fit$h)$sigma
		FPF.aux <- matrix(1 - ecdf(h.residuals)(aux), nrow = length(TPF))
	}
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
