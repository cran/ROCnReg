summary.cROC <-
function(object, ...) {
	res <- list()
	res$call <- object$call
	method <- switch(class(object)[1], "cROC.kernel" = "Conditional ROC curve - Kernel-based", 
								  	   "cROC.bnp" = "Conditional ROC curve - Bayesian nonparametric", 
								   	   "cROC.sp" = "Conditional ROC curve - semiparametric")

	res$method <- method

	if(class(object)[1] == "cROC.kernel") {
		m <- matrix(ncol = 2, nrow = 1, dimnames = list(c("Bandwidth:"), c("Group H", "Group D")))
		m[1,] <- c(sprintf("%.6f", object$fit$h$bw.mean$bw), sprintf("%.6f", object$fit$d$bw.mean$bw))
		res$kernel.regfun$bw <- m
		attr(res$kernel.regfun, "pregtype") <- paste0("\nKernel Estimator: ", object$fit$h$bw.mean$pregtype)
		attr(res$kernel.regfun, "pmethod") <- paste0("\nBandwidth Selection Method: ", object$fit$h$bw.mean$pmethod)
		attr(res$kernel.regfun, "pckertype") <- paste0("\nContinuous Kernel Type: ", object$fit$h$bw.mean$klist$x$pckertype)

		m <- matrix(ncol = 2, nrow = 1, dimnames = list(c("Bandwidth:"), c("Group H", "Group D")))
		m[1,] <- c(sprintf("%.6f", object$fit$h$bw.var$bw), sprintf("%.6f", object$fit$d$bw.var$bw))
		
		res$kernel.varfun$bw <- m
		attr(res$kernel.varfun, "pregtype") <- paste0("\nKernel Estimator: ", object$fit$h$bw.var$pregtype)
		attr(res$kernel.varfun, "pmethod") <- paste0("\nBandwidth Selection Method: ", object$fit$h$bw.var$pmethod)
		attr(res$kernel.varfun, "pckertype") <- paste0("\nContinuous Kernel Type: ", object$fit$h$bw.var$klist$x$pckertype)
	}

	if(class(object)[1] == "cROC.sp") {
		if(ncol(object$coeff$h) == 3) {
			colnames(object$coeff$h) <- c("Estimate", "Quantile 2.5%", "Quantile 97.5%")
		} else {
			colnames(object$coeff$h) <- c("Estimate")
		}
		res$sp.coeff$h <- object$coeff$h
		if(ncol(object$coeff$d) == 3) {
			colnames(object$coeff$d) <- c("Estimate", "Quantile 2.5%", "Quantile 97.5%")
		} else {
			colnames(object$coeff$d) <- c("Estimate")
		}
		res$sp.coeff$d <- object$coeff$d

		if(ncol(object$coeff$ROC) == 3) {
			colnames(object$coeff$ROC) <- c("Estimate", "Quantile 2.5%", "Quantile 97.5%")
		} else {
			colnames(object$coeff$ROC) <- c("Estimate")
		}
		res$sp.coeff$ROC <- object$coeff$ROC

		col.names <- c("Group H", "Group D")
		row.names <- c("AIC", "BIC")
		m <- matrix(ncol = length(col.names), nrow = length(row.names), dimnames = list(row.names, col.names))
		n0 <- nrow((object$data[object$data[,object$group] == object$tag.h,])[!object$missing.ind$h,])
		n1 <- nrow((object$data[object$data[,object$group] != object$tag.h,])[!object$missing.ind$d,])
		m[1,] <- c(sprintf("%.3f", AIC(object$fit$h)), sprintf("%.3f", AIC(object$fit$d)))
		m[2,] <- c(sprintf("%.3f", AIC(object$fit$h, k = log(n0))), sprintf("%.3f", AIC(object$fit$d, k = log(n1))))

		res$sp.msc <- m
	}
	if(class(object)[1] == "cROC.bnp" && object$prior$h$L == 1 && object$prior$d$L == 1) {
		beta.h <- object$fit$h$beta[,object$fit$h$mm$paracoeff, drop = FALSE]
		beta.d <- object$fit$d$beta[,object$fit$d$mm$paracoeff, drop = FALSE]

		m <- matrix(ncol = 3, nrow = ncol(beta.h), dimnames = list(colnames(beta.h), c("Post. mean", "Post. quantile 2.5%", "Post. quantile 97.5%")))
		for(i in 1:ncol(beta.h)) {
			#m[i,] <- c(sprintf("%.4f", mean(beta.h[,i], na.rm = TRUE)), sprintf("%.4f", quantile(beta.h[,i], 0.025, na.rm = TRUE)), sprintf("%.4f", quantile(beta.h[,i], 0.975, na.rm = TRUE)))
			m[i,] <- c(mean(beta.h[,i], na.rm = TRUE), quantile(beta.h[,i], 0.025, na.rm = TRUE), quantile(beta.h[,i], 0.975, na.rm = TRUE))
		}
		res$bnp.coeff$h <- m

		m <- matrix(ncol = 3, nrow = ncol(beta.d), dimnames = list(colnames(beta.d), c("Post. mean", "Post. quantile 2.5%", "Post. quantile 97.5%")))
		for(i in 1:ncol(beta.d)) {
			#m[i,] <- c(sprintf("%.4f", mean(beta.d[,i], na.rm = TRUE)), sprintf("%.4f", quantile(beta.d[,i], 0.025, na.rm = TRUE)), sprintf("%.4f", quantile(beta.d[,i], 0.975, na.rm = TRUE)))
			m[i,] <- c(mean(beta.d[,i], na.rm = TRUE), quantile(beta.d[,i], 0.025, na.rm = TRUE), quantile(beta.d[,i], 0.975, na.rm = TRUE))
		}
		res$bnp.coeff$d <- m
		
		coeff.h <- colnames(beta.h)
		coeff.d <- colnames(beta.d)

		coeffs <- c(coeff.h, coeff.d[is.na(match(coeff.d, coeff.h))])
		
		beta.h.ROC <- beta.d.ROC <- matrix(0, ncol = length(coeffs), nrow = nrow(beta.h))
		colnames(beta.h.ROC) <- colnames(beta.d.ROC) <- coeffs
		
		beta.h.ROC[,match(coeff.h, coeffs)] <- beta.h
		beta.d.ROC[,match(coeff.d, coeffs)] <- beta.d

		beta.ROC <- cbind((beta.h.ROC - beta.d.ROC)/object$fit$d$sd, "b" = object$fit$h$sd/object$fit$d$sd)
		m <- matrix(ncol = 3, nrow = ncol(beta.ROC), dimnames = list(colnames(beta.ROC), c("Post. mean", "Post. quantile 2.5%", "Post. quantile 97.5%")))
		for(i in 1:ncol(beta.ROC)) {
			#m[i,] <- c(sprintf("%.4f", mean(beta.ROC[,i], na.rm = TRUE)), sprintf("%.4f", quantile(beta.ROC[,i], 0.025, na.rm = TRUE)), sprintf("%.4f", quantile(beta.ROC[,i], 0.975, na.rm = TRUE)))
			m[i,] <- c(mean(beta.ROC[,i], na.rm = TRUE), quantile(beta.ROC[,i], 0.025, na.rm = TRUE), quantile(beta.ROC[,i], 0.975, na.rm = TRUE))
		}
		res$bnp.coeff$ROC <- m
	}

	waic <- !is.null(object$WAIC)
	lpml <- !is.null(object$lpml)
	dic  <- !is.null(object$DIC)

	if(waic | lpml | dic) {
		col.names <- c("Group H", "Group D")
		row.names <- NULL
		m <- matrix(ncol = length(col.names), nrow = ifelse(waic, 2, 0) + ifelse(lpml, 1, 0) + ifelse(dic, 2, 0))
		i <- 1
		if(waic) {
			row.names <- c(row.names, "WAIC", "WAIC (Penalty)")
			m[i,] <-  c(sprintf("%.3f", object$WAIC$h$WAIC), sprintf("%.3f", object$WAIC$d$WAIC)) 
			m[i+1,] <- c(sprintf("%.3f", object$WAIC$h$pW), sprintf("%.3f", object$WAIC$d$pW))
			i <- i + 2
		}
		if(lpml) {
			row.names <- c(row.names, "LPML")
			m[i,]   <-  c(sprintf("%.3f", object$lpml$h$lpml), sprintf("%.3f", object$lpml$d$lpml))
			i <- i + 1
		}
		if(dic) {
			row.names <- c(row.names, "DIC", "DIC (Penalty)")
			m[i,] <-  c(sprintf("%.3f", object$DIC$h$DIC), sprintf("%.3f", object$DIC$d$DIC)) 
			m[i+1,] <- c(sprintf("%.3f", object$DIC$h$pD), sprintf("%.3f", object$DIC$d$pD))
		}
		colnames(m) <- col.names
		rownames(m) <- row.names

		res$bmsc <- m
	}
	m <- matrix(ncol = 2, nrow = 2, dimnames = list(c("Number of observations", "Number of missing data"), c("Group H", "Group D")))
	m[1,] <- c(sprintf("%.0f", nrow(object$data[object$data[,object$group] == object$tag.h,])), sprintf("%.0f", nrow(object$data[object$data[,object$group] != object$tag.h,])))
	m[2,] <- c(sprintf("%.0f", sum(object$missing.ind$h)), sprintf("%.0f", sum(object$missing.ind$d)))
	res$sz <- m

	print.summary.cROC(res, ...)
	invisible(res)   	   	
}
