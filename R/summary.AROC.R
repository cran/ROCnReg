summary.AROC <-
function(object, ...) {
	res <- list()
	res$call <- object$call
	
	method <- switch(class(object)[1], "AROC.kernel" = "AROC Kernel-based", 
									   "AROC.bnp" = "AROC Bayesian nonparametric", 
									   "AROC.sp" = "AROC semiparametric")
	res$method <- method
	auc_aauc <- "Area under the covariate-adjusted ROC curve"
	if(length(object$AUC) == 3) {
		AUC <- paste0(auc_aauc, ": ", paste(round(object$AUC[1], 3), " (", round(object$AUC[2], 3),"",", ", round(object$AUC[3], 3),")", sep = ""))
	} else {
		AUC <- paste0(auc_aauc, ": ", round(object$AUC[1], 3))
	}
	res$AUC <- AUC

	if(!is.null(object$pAUC)) {
		p_auc_aauc <- ifelse(attr(object$pAUC, "focus") == "FPF", "Partial area under the covariate-adjusted ROC curve", "Partial area under the specificity covariate-adjusted ROC curve")
		p_auc_aauc <- paste0(p_auc_aauc, ifelse(attr(object$pAUC, "focus") == "FPF", " (FPF = ", " (Se = "), attr(object$pAUC, "value"), ")")

		if(length(object$pAUC) == 3) {
			pAUC <- paste0(p_auc_aauc, ": ", paste(round(object$pAUC[1], 3), " (", round(object$pAUC[2], 3),"",", ", round(object$pAUC[3], 3),")", sep = ""))
		} else {
			pAUC <- paste0(p_auc_aauc, ": ", round(object$pAUC[1], 3))
		}
		res$pAUC <- pAUC
	}

	if(class(object)[1] == "AROC.kernel") {
		m <- matrix(ncol = 1, nrow = 1, dimnames = list(c("Bandwidth:"), c("Group H")))
		m[1,] <- sprintf("%.6f", object$fit$bw.mean$bw)
		res$kernel.regfun$bw <- m
		attr(res$kernel.regfun, "pregtype") <- paste0("\nKernel Estimator: ", object$fit$bw.mean$pregtype)
		attr(res$kernel.regfun, "pmethod") <- paste0("\nBandwidth Selection Method: ", object$fit$bw.mean$pmethod)
		attr(res$kernel.regfun, "pckertype") <- paste0("\nContinuous Kernel Type: ", object$fit$bw.mean$klist$x$pckertype)
		
		m <- matrix(ncol = 1, nrow = 1, dimnames = list(c("Bandwidth:"), c("Group H")))
		m[1,] <- sprintf("%.6f", object$fit$bw.var$bw)
		res$kernel.varfun$bw <- m
		attr(res$kernel.varfun, "pregtype") <- paste0("\nKernel Estimator: ", object$fit$bw.var$pregtype)
		attr(res$kernel.varfun, "pmethod") <- paste0("\nBandwidth Selection Method: ", object$fit$bw.var$pmethod)
		attr(res$kernel.varfun, "pckertype") <- paste0("\nContinuous Kernel Type: ", object$fit$bw.var$klist$x$pckertype)
	}

	if(class(object)[1] == "AROC.sp") {
		if(ncol(object$coeff) == 3) {
			colnames(object$coeff) <- c("Estimate", "Quantile 2.5%", "Quantile 97.5%")
		} else {
			colnames(object$coeff) <- c("Estimate")
		}
		res$sp.coeff <- object$coeff
		col.names <- c("Group H")
		row.names <- c("AIC", "BIC")
		m <- matrix(ncol = length(col.names), nrow = length(row.names), dimnames = list(row.names, col.names))
		n0 <- nrow((object$data[object$data[,object$group] == object$tag.h,])[!object$missing.ind$h,])
		m[1,] <- sprintf("%.3f", AIC(object$fit))
		m[2,] <- sprintf("%.3f", AIC(object$fit, k = log(n0)))
		res$sp.msc <- m
	}
	if(class(object)[1] == "AROC.bnp" &&  object$prior$L == 1) {
		beta.h <- object$fit$beta[,object$fit$mm$paracoeff,drop = FALSE]
		m <- matrix(ncol = 3, nrow = ncol(beta.h), dimnames = list(colnames(beta.h), c("Post. mean", "Post. quantile 2.5%", "Post. quantile 97.5%")))
		for(i in 1:ncol(beta.h)) {
			#m[i,] <- c(sprintf("%.4f", mean(beta.h[,i], na.rm = TRUE)), sprintf("%.4f", quantile(beta.h[,i], 0.025, na.rm = TRUE)), sprintf("%.4f", quantile(beta.h[,i], 0.975, na.rm = TRUE)))
			m[i,] <- c(mean(beta.h[,i], na.rm = TRUE), quantile(beta.h[,i], 0.025, na.rm = TRUE), quantile(beta.h[,i], 0.975, na.rm = TRUE))
		}
		res$bnp.coeff <- m
	}

	waic <- !is.null(object$WAIC)
	lpml <- !is.null(object$lpml)
	dic  <- !is.null(object$DIC)

	if(waic | lpml | dic) {
		col.names <- c("Group H")
		row.names <- NULL
		m <- matrix(ncol = length(col.names), nrow = ifelse(waic, 2, 0) + ifelse(lpml, 1, 0) + ifelse(dic, 2, 0))
		i <- 1
		if(waic) {
			row.names <- c(row.names, "WAIC", "WAIC (Penalty)")
			m[i,1] <-  sprintf("%.3f", object$WAIC$WAIC)
			m[i+1,1] <- sprintf("%.3f", object$WAIC$pW)
			i <- i+2
		}
		if(lpml) {
			row.names <- c(row.names, "LPML")
			m[i,1]   <-  sprintf("%.3f", object$lpml$lpml)
			i <- i + 1
		}
		if(dic) {
			row.names <- c(row.names, "DIC", "DIC (Penalty)")
			m[i,1]   <-  sprintf("%.3f", object$DIC$DIC)
			m[i+1,1] <- sprintf("%.3f", object$DIC$pD)
		}
		colnames(m) <- col.names
		rownames(m) <- row.names
		res$bmsc <- m
	}
	m <- matrix(ncol = 2, nrow = 2, dimnames = list(c("Number of observations", "Number of missing data"), c("Group H", "Group D")))
	m[1,] <- c(sprintf("%.0f", nrow(object$data[object$data[,object$group] == object$tag.h,])), sprintf("%.0f", nrow(object$data[object$data[,object$group] != object$tag.h,])))
	m[2,] <- c(sprintf("%.0f", sum(object$missing.ind$h)), sprintf("%.0f", sum(object$missing.ind$d)))
	res$sz <- m
	
	print.summary.AROC(res, ...)
	invisible(res)	   	   	
}
