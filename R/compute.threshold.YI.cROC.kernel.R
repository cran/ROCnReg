compute.threshold.YI.cROC.kernel <-
function(object, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	if(class(object)[1] != "cROC.kernel") {
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

	# Compute F_D|X and F_{\bar{D}}|X
	fit1.mean.new <- npreg(object$fit$d$bw.mean, exdat = xp, residuals = TRUE)
	fit1.var.new <- npreg(object$fit$d$bw.var, exdat = xp, residuals = TRUE)
	res1p <- object$fit$d$fit.mean$resid/sqrt(object$fit$d$fit.var$mean)

	fit0.mean.new <- npreg(object$fit$h$bw.mean, exdat = xp, residuals = TRUE)
	fit0.var.new <- npreg(object$fit$h$bw.var, exdat = xp, residuals = TRUE)
	res0p <- object$fit$h$fit.mean$resid/sqrt(object$fit$h$fit.var$mean)

	npred <- length(xp)

	y0 <- (object$data[object$data[,object$group] == object$tag.h,])[!object$missing.ind$h, object$marker]
    y1 <- (object$data[object$data[,object$group] != object$tag.h,])[!object$missing.ind$d, object$marker]
  	
  	n0 <- length(y0)
  	n1 <- length(y1)

  	grid  <- seq(min(c(y0, y1), na.rm = TRUE) - 1, max(c(y0, y1), na.rm = TRUE) + 1, length = max(500, c(n0,n1)))
	
	#grid  <- seq(min(object$data[, object$marker], na.rm = TRUE) - 1, max(object$data[, object$marker], na.rm = TRUE) + 1, length = 500)
	ngrid <- length(grid)

	F0 <- F1 <- matrix(0, nrow = ngrid, ncol = npred)
	thresholds.s <- YI.s <- TPF.s <- FPF.s <- vector(length = npred)

	for(l in 1:npred) {
		F0[,l] <- ecdf(res0p)((grid - fit0.mean.new$mean[l])/sqrt(fit0.var.new$mean[l]))
		F1[,l] <- ecdf(res1p)((grid - fit1.mean.new$mean[l])/sqrt(fit1.var.new$mean[l]))

		difbb <- abs(F0[,l] -  F1[,l])
		thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])

		YI.s[l] <- max(difbb)
		TPF.s[l] <- 1 - ecdf(res1p)((thresholds.s[l] - fit1.mean.new$mean[l])/sqrt(fit1.var.new$mean[l]))
		FPF.s[l] <- 1 - ecdf(res0p)((thresholds.s[l] - fit0.mean.new$mean[l])/sqrt(fit0.var.new$mean[l]))
	}
	res <- list()
	res$call <- match.call()
	res$newdata <- newdata
	res$thresholds <- thresholds.s
	res$YI  <- YI.s
	res$FPF <- FPF.s
	res$TPF <- TPF.s
	res
  
}
