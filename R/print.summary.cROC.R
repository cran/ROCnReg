print.summary.cROC <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", x$method))
	cat("\n----------------------------------------------------------")

	if(!is.null(x$kernel.regfun)) {
		cat("\n\nRegression functions:\n\n")
		print(x$kernel.regfun$bw, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
		cat(attr(x$kernel.regfun, "pregtype"))
		cat(attr(x$kernel.regfun, "pmethod"))
		cat(attr(x$kernel.regfun, "pckertype"))
	}
	if(!is.null(x$kernel.varfun)) {
		cat("\n\nVariance functions:\n\n")
		print(x$kernel.varfun$bw, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
		cat(attr(x$kernel.varfun, "pregtype"))
		cat(attr(x$kernel.varfun, "pmethod"))
		cat(attr(x$kernel.varfun, "pckertype"))
	}
	if(!is.null(x$sp.coeff)) {
		cat("\n\nParametric coefficients")
		cat("\nGroup H:\n")
		#print(x$sp.coeff$h, quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
		print(format(round(x$sp.coeff$h, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
		cat("\n\nGroup D:\n")
		#print(x$sp.coeff$d, quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
		print(format(round(x$sp.coeff$d, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
		cat("\n\nROC curve:\n")
		#print(x$sp.coeff$ROC, quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
		print(format(round(x$sp.coeff$ROC, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
	}
	if(!is.null(x$sp.msc)) {
		cat("\n\nModel selection criteria:\n")
		print(x$sp.msc, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}

	if(!is.null(x$bnp.coeff)) {
		cat("\n\nParametric coefficients")
		cat("\nGroup H:\n")
		print(format(round(x$bnp.coeff$h, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)

		cat("\n\nGroup D:\n")
		print(format(round(x$bnp.coeff$d, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)

		cat("\n\nROC curve:\n")
		print(format(round(x$bnp.coeff$ROC, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
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
