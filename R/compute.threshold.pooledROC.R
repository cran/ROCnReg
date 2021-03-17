compute.threshold.pooledROC <-
function(object, criterion = c("FPF", "YI"), FPF, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	object.class <- class(object)[1]
	if (!(object.class %in% c("pooledROC.BB", "pooledROC.emp", "pooledROC.emp", "pooledROC.kernel", "pooledROC.dpm"))) {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}

	criterion <- match.arg(criterion)
	if(criterion == "FPF" & missing(FPF)) {
		stop(paste0("The vector of FPF at which to calculate the threshold values should be specified"))	
	}
	if(ci.level <= 0 || ci.level >= 1) {
        stop("The ci.level should be between 0 and 1")
    }
	method <- paste0("compute.threshold.", criterion, ".", object.class)
	if(criterion == "YI") {
		res <- eval(parse(text = method))(object = object, ci.level = ci.level, parallel = parallel, ncpus = ncpus, cl = cl)
	} else {
		res <- eval(parse(text = method))(object = object, FPF = FPF, ci.level = ci.level, parallel = parallel, ncpus = ncpus, cl = cl)
	}
	res$call <- match.call()
	res
}
