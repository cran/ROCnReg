print.AROC <-
function(x, ...) {
	method <- switch(class(x)[1], "AROC.kernel" = "AROC Kernel-based", "AROC.bnp" = "AROC Bayesian nonparametric", "AROC.bsp" = "AROC Bayesian semiparametric", "AROC.sp" = "AROC semiparametric")

	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", method))
	cat("\n----------------------------------------------\n")

	auc_aauc <- "Area under the covariate-adjusted ROC curve"
	if(length(x$AUC) == 3) {
		AUC <- paste0(auc_aauc, ": ", paste(round(x$AUC[1], 3), " (", round(x$AUC[2], 3),"",", ", round(x$AUC[3], 3),")", sep = ""))
	} else {
		AUC <- paste0(auc_aauc, ": ", round(x$AUC[1], 3))
	}
	cat(AUC, "\n")

	if(!is.null(x$pAUC)) {
		#p_auc_aauc <- "Partial area under the covariate-adjusted ROC curve"
		#p_auc_aauc <- paste0(p_auc_aauc, " (FPF = ", attr(x$pAUC, "value"), ")")

		p_auc_aauc <- ifelse(attr(x$pAUC, "focus") == "FPF", "Partial area under the covariate-adjusted ROC curve", "Partial area under the specificity covariate-adjusted ROC curve")
		p_auc_aauc <- paste0(p_auc_aauc, ifelse(attr(x$pAUC, "focus") == "FPF", " (FPF = ", " (Se = "), attr(x$pAUC, "value"), ")")


		if(length(x$pAUC) == 3) {
			pAUC <- paste0(p_auc_aauc, ": ", paste(round(x$pAUC[1], 3), " (", round(x$pAUC[2], 3),"",", ", round(x$pAUC[3], 3),")", sep = ""))
		} else {
			pAUC <- paste0(p_auc_aauc, ": ", round(x$pAUC[1], 3))
		}
		cat(paste0(pAUC, "\n"))
	}
	invisible(x)
}
