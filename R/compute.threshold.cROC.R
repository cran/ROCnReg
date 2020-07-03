compute.threshold.cROC <-
function(object, criterion = c("FPF", "YI"), FPF, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	object.class <- class(object)[1]
	if (!(object.class %in% c("cROC.bnp", "cROC.kernel", "cROC.sp"))) {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}

	criterion <- match.arg(criterion)
	if(criterion == "FPF" & missing(FPF)) {
		stop(paste0("The vector of FPF at which to calculate the threshold values should be specified"))	
	}
	method <- paste0("compute.threshold.", criterion, ".", object.class)
	if(criterion == "YI") {
		res <- eval(parse(text = method))(object = object, newdata = newdata, parallel = parallel, ncpus = ncpus, cl = cl)
	} else {
		res <- eval(parse(text = method))(object = object, FPF = FPF, newdata = newdata, parallel = parallel, ncpus = ncpus, cl = cl)
	}
	res$call <- match.call()
	res
}
