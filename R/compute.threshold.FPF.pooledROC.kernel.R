compute.threshold.FPF.pooledROC.kernel <-
function(object, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "pooledROC.kernel") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}
	thresholds <- TPF <- vector(length = length(FPF))
	for(i in 1:length(FPF)) {
		thresholds[i] <- qFk(1-FPF[i], y = object$marker$h[!object$missing.ind$h], h = object$bws$h)
		TPF[i] <- 1-Gk(thresholds[i], y = object$marker$d[!object$missing.ind$d], h = object$bws$d)
	}
	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res
}
