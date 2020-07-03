compute.threshold.FPF.pooledROC.emp <-
function(object, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "pooledROC.emp") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}
	
	F1emp <- ecdf(object$marker$d[!object$missing.ind$d])
	thresholds <- quantile(object$marker$h[!object$missing.ind$h], 1 - FPF, type = 1)
	TPF <- 1 - F1emp(thresholds)
	
	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res

}
