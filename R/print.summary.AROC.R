print.summary.AROC <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", x$method))
	cat("\n----------------------------------------------\n")
	if(!is.null(x$AUC)) {
		cat(x$AUC)
	}
	if(!is.null(x$pAUC)) {
		cat(paste0("\n", x$pAUC))
	}
	if(!is.null(x$kernel.regfun)) {
		cat("\n\nRegression function:\n\n")
		print(x$kernel.regfun$bw, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
		cat(attr(x$kernel.regfun, "pregtype"))
		cat(attr(x$kernel.regfun, "pmethod"))
		cat(attr(x$kernel.regfun, "pckertype"))
	}
	if(!is.null(x$kernel.varfun)) {
		cat("\n\nVariance function:\n\n")
		print(x$kernel.varfun$bw, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
		aux <- attr(x$kernel.varfun, "pregtype")
		cat(aux)
		cat(attr(x$kernel.varfun, "pmethod"))
		cat(attr(x$kernel.varfun, "pckertype"))
	}
	if(!is.null(x$sp.coeff)) {
		cat("\n\nParametric coefficients (Group H):\n")
		#print(x$sp.coeff, quote = FALSE, right = TRUE, na.print = "", print.gap = 5, digits = digits)
		print(format(round(x$sp.coeff, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 4, digits = digits)
	}
	if(!is.null(x$sp.msc)) {
		cat("\n\nModel selection criteria:\n")
		print(x$sp.msc, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	if(!is.null(x$bnp.coeff)) {
		cat("\n\nParametric coefficients (Group H):\n")
		print(format(round(x$bnp.coeff, digits), digits = digits), quote = FALSE, right = TRUE, na.print = "", print.gap = 5, digits = digits)
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
