print.cROC <-
function(x, ...) {
	method <- switch(class(x)[1], "cROC.kernel" = "Conditional ROC curve - Kernel-based", 
								  "cROC.bnp" = "Conditional ROC curve - Bayesian nonparametric", 
								  "cROC.bsp" = "Conditional ROC curve - Bayesian semiparametric", 
								  "cROC.sp" = "Conditional ROC curve - semiparametric")

	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", method))
	cat("\n----------------------------------------------------------\n")

	invisible(x)   
}
