compute.threshold.TPF.pooledROC.kernel <-
function(object, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "pooledROC.kernel") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}
	thresholds <- FPF <- vector(length = length(FPF))
	for(i in 1:length(TPF)) {
		thresholds[i] <- qFk(1-TPF[i], y = object$marker$d[!object$missing.ind$d], h = object$bws$d)
		FPF[i] <- 1-Gk(thresholds[i], y = object$marker$h[!object$missing.ind$h], h = object$bws$h)
	}
	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res
}
