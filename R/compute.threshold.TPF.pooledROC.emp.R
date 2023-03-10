compute.threshold.TPF.pooledROC.emp <-
function(object, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "pooledROC.emp") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}
	
	F0emp <- ecdf(object$marker$h[!object$missing.ind$h])
	thresholds <- quantile(object$marker$d[!object$missing.ind$d], 1 - TPF, type = 1)
	FPF <- 1 - F0emp(thresholds)
	
	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res

}
