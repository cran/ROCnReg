print.summary.pooledROC <-
function(x,...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", x$method))
	cat("\n----------------------------------------------\n")
	if(!is.null(x$AUC)) {
		cat(x$AUC)
	}
	if(!is.null(x$pAUC)) {
		cat(paste0("\n", x$pAUC))
	}
	if(!is.null(x$bws)) {
		cat("\n\n")
		print(x$bws, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	if(!is.null(x$bw)) {
		cat("\n\n")
		cat(x$bw)
	}
	if(!is.null(x$bmsc)) {
		cat("\n\nModel selection criteria:\n")
		print(x$bmsc, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	if(!is.null(x$sz)) {
		cat("\n\nSample sizes:\n")
		print(x$sz, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	invisible(x)
}
