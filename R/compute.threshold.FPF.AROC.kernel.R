compute.threshold.FPF.AROC.kernel <-
function(object, newdata, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "AROC.kernel") {
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

	res.aux <- compute.threshold.FPF.kernel(object = object$fit, newdata = xp, FPF = FPF)

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
	res$newdata <- newdata
	res

}
