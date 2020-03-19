compute.threshold.pooledROC <-
function(object, criterion = c("FPF", "YI"), FPF) {
	object.class <- class(object)[2]
	if (!(object.class %in% c("pooledROC.BB", "pooledROC.emp", "pooledROC.emp", "pooledROC.kernel", "pooledROC.dpm"))) {
		stop(paste0("This function can not be used for this object class: ", class(object)[2]))
	}

	criterion <- match.arg(criterion)
	if(criterion == "FPF" & missing(FPF)) {
		stop(paste0("The vector of FPF at which to calculate the threshold values should be specified"))	
	}
	method <- paste0("compute.threshold.", criterion, ".", object.class)
	if(criterion == "YI") {
		res <- eval(parse(text = method))(object = object)
	} else {
		res <- eval(parse(text = method))(object = object, FPF = FPF)
	}
	res$call <- match.call()
	res
}