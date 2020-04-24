summary.pooledROC <-
function(object, ...) {
	res <- list()	
	res$call <- object$call
	method <- switch(class(object)[1], "pooledROC.BB" = "Pooled ROC curve - Bayesian bootstrap", "pooledROC.emp" = "Pooled ROC curve - Empirical", "pooledROC.emp" = "Pooled ROC curve - Empirical", "pooledROC.kernel" = "Pooled ROC curve - Kernel-based", "pooledROC.dpm" = "Pooled ROC curve - Bayesian DPM")
	res$method <- method
	# AUC and pAUC
	auc_aauc <- "Area under the pooled ROC curve"
	if(length(object$AUC) == 3) {
		legend.text <- paste0(auc_aauc, ": ", paste(round(object$AUC[1], 3), " (", round(object$AUC[2], 3),"",", ", round(object$AUC[3], 3),")", sep = ""))
	} else {
		legend.text <- paste0(auc_aauc, ": ", round(object$AUC[1], 3))
	}
	res$AUC <- legend.text
	if(!is.null(object$pAUC)) {
		p_auc_aauc <- ifelse(attr(object$pAUC, "focus") == "FPF", "Partial area under the pooled ROC curve", "Partial area under the specificity pooled ROC curve")
		p_auc_aauc <- paste0(p_auc_aauc, ifelse(attr(object$pAUC, "focus") == "FPF", " (FPF = ", " (Se = "), attr(object$pAUC, "value"), ")")

		if(length(object$pAUC) == 3) {
			legend.text <- paste0(p_auc_aauc, ": ", paste(round(object$pAUC[1], 3), " (", round(object$pAUC[2], 3),"",", ", round(object$pAUC[3], 3),")", sep = ""))
		} else {
			legend.text <- paste0(p_auc_aauc, ": ", round(object$pAUC[1], 3))
		}
		res$pAUC <- legend.text
	}

	if(class(object)[1] == "pooledROC.kernel") {
		m <- matrix(ncol = 2, nrow = 1, dimnames = list(c("Bandwidths:"), c("Group H", "Group D")))
		m[1,] <- c(sprintf("%.3f", object$bws$h), sprintf("%.3f", object$bws$d))
		res$bws <- m
		res$bw <- paste0("\nBandwidth Selection Method: ", switch(object$bw, "SRT" = "Silverman's rule-of-thumb", "UCV" = "Unbiased cross-validation"))
	}

	waic <- any(class(object) %in% "pooledROC.dpm") & !is.null(object$WAIC)
	lpml <- any(class(object) %in% "pooledROC.dpm") & !is.null(object$lpml)
	dic  <- any(class(object) %in% "pooledROC.dpm") & !is.null(object$DIC)

	if(waic | lpml | dic) {
		col.names <- c("Group H", "Group D")
		row.names <- NULL
		m <- matrix(ncol = length(col.names), nrow = ifelse(waic, 2, 0) + ifelse(lpml, 1, 0) + ifelse(dic, 2, 0))
		i <- 1
		if(waic) {
			row.names <- c(row.names, "WAIC", "WAIC (Penalty)")
			m[i,] <-  c(sprintf("%.3f", object$WAIC$h$WAIC), sprintf("%.3f", object$WAIC$d$WAIC)) 
			m[i+1,] <- c(sprintf("%.3f", object$WAIC$h$pW), sprintf("%.3f", object$WAIC$d$pW))
			i <- i+2
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
	m[1,] <- c(sprintf("%.0f", length(object$marker$h)), sprintf("%.0f", length(object$marker$d)))
	m[2,] <- c(sprintf("%.0f", sum(object$missing.ind$h)), sprintf("%.0f", sum(object$missing.ind$d)))
	res$sz <- m
	print.summary.pooledROC(res)
	invisible(res)	   	
}
