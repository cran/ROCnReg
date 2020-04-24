print.pooledROC <-
function(x, ...) {
	method <- switch(class(x)[1], "pooledROC.BB" = "Pooled ROC curve - Bayesian bootstrap", "pooledROC.emp" = "Pooled ROC curve - Empirical", "pooledROC.emp" = "Pooled ROC curve - Empirical", "pooledROC.kernel" = "Pooled ROC curve - Kernel-based", "pooledROC.dpm" = "Pooled ROC curve - Bayesian DPM")
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", method))
	cat("\n----------------------------------------------\n")
	# AUC and pAUC
	auc_aauc <- "Area under the pooled ROC curve"
	if(length(x$AUC) == 3) {
		legend.text <- paste0(auc_aauc, ": ", paste(round(x$AUC[1], 3), " (", round(x$AUC[2], 3),"",", ", round(x$AUC[3], 3),")", sep = ""))
	} else {
		legend.text <- paste0(auc_aauc, ": ", round(x$AUC[1], 3))
	}
	cat(legend.text, "\n")

	if(!is.null(x$pAUC)) {
		p_auc_aauc <- ifelse(attr(x$pAUC, "focus") == "FPF", "Partial area under the pooled ROC curve", "Partial area under the specificity pooled ROC curve")
		p_auc_aauc <- paste0(p_auc_aauc, ifelse(attr(x$pAUC, "focus") == "FPF", " (FPF = ", " (Se = "), attr(x$pAUC, "value"), ")")

		if(length(x$pAUC) == 3) {
			legend.text <- paste0(p_auc_aauc, ": ", paste(round(x$pAUC[1], 3), " (", round(x$pAUC[2], 3),"",", ", round(x$pAUC[3], 3),")", sep = ""))
		} else {
			legend.text <- paste0(p_auc_aauc, ": ", round(x$pAUC[1], 3))
		}
		cat(paste0(legend.text, "\n"))
	}
	invisible(x)
}
