compute.threshold.FPF.cROC.sp <-
function(object, newdata, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
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
		newdata <- na.omit(newdata[,names.cov,drop = FALSE])
	}  
	res.aux <- compute.threshold.FPF.sp(object = list(est.cdf = object$est.cdf, fit = object$fit$h), newdata = newdata, FPF = FPF)
	
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
	
	fit.new <- predict(object$fit$d, newdata = newdata)
	aux <- t(t(res.aux$thresholds) - fit.new)
	aux <- t(t(aux)/summary(object$fit$d)$sigma)

	if(object$est.cdf == "normal") {
		TPF.aux <- 1 - pnorm(aux)
	} else {
		d.residuals <- object$fit$d$residuals/summary(object$fit$d)$sigma
		TPF.aux <- matrix(1 - ecdf(d.residuals)(aux), nrow = length(FPF))
	}
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
