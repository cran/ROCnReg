compute.threshold.YI.AROC.kernel <-
function(object, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "AROC.kernel") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}
  # Newdata
	names.cov <- object$covariate

	if(!missing(newdata) && !inherits(newdata, "data.frame"))
		stop("Newdata must be a data frame")
	if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
		stop("Not all needed variables are supplied in newdata") 

	if(missing(newdata)) {
		newdata <- cROCData(object$data, names.cov, object$group)
	} else {
		newdata <- na.omit(newdata[,names.cov,drop = FALSE])
	}

	xp <- newdata[,names.cov]

	p <- seq(0, 1, length = 500)
	np <- length(p)
	npred <- length(xp)

	# Compute the AROC
	data.d <- (object$data[object$data[,object$group] != object$tag.h,])[!object$missing.ind$d,]
	x1 <- data.d[,object$covariate]
	y1 <- data.d[,object$marker]
	n1 <- length(y1)

	# Evaluate the model in the diseased population, and compute the AROC
	fit.mean.d.p <- npreg(object$fit$bw.mean, exdat = x1,residuals = TRUE)
	fit.var.d.p <- npreg(object$fit$bw.var, exdat = x1, residuals = TRUE)

	res0p <- object$fit$fit.mean$resid/sqrt(object$fit$fit.var$mean)
	F0res <- ecdf(res0p)

	u1 <- 1 - F0res((y1 - fit.mean.d.p$mean)/sqrt(fit.var.d.p$mean))
	AROC <- numeric(np)
	for(i in 1:np){
		AROC[i] <- sum(u1 <= p[i])/n1
	}
	# Compute YI and associated threshold values
	difbb <- AROC -  p
	FPF <- mean(p[which(difbb == max(difbb))])
	YI <- max(difbb)

	fit0.mean.new <- npreg(object$fit$bw.mean, exdat = xp, residuals = TRUE)
	fit0.var.new <- npreg(object$fit$bw.var, exdat = xp, residuals = TRUE)
	thresholds <- fit0.mean.new$mean + sqrt(fit0.var.new$mean)*quantile(res0p, 1 - FPF, type = 1)

	res <- list()
	res$call <- match.call()
	res$newdata <- newdata
	res$thresholds <- thresholds
	res$YI  <- YI
	res$FPF <- FPF
	res
  
}
