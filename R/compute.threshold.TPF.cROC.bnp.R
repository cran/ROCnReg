compute.threshold.TPF.cROC.bnp <-
function(object, newdata, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    if(class(object)[1] != "cROC.bnp") {
        stop(paste0("This function cannot be used for this object class: ", class(object)[1]))
    }
    # Newdata
    #names.cov.h <- all.vars(object$fit$h$formula)[-1]
    #names.cov.d <- all.vars(object$fit$d$formula)[-1]
    #names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])

    names.cov.h <- get_vars_formula(object$fit$h$formula)
    names.cov.d <- get_vars_formula(object$fit$d$formula)
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
    
    thresholds <- compute.threshold.TPF.bnp(object_h = object$fit$h, object_d = object$fit$d, newdata = newdata, TPF = TPF, ci.level = ci.level, parallel = parallel, ncpus = ncpus, cl = cl)
        
    res <- list()
    res$newdata <- newdata
    res$thresholds <- thresholds$thresholds
    res$FPF <- thresholds$FPF
    res$TPF <- TPF
    res
}
