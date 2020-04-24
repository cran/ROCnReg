plot.AROC <-
function(x, main = NULL, ...) {
    dots <- list(...)

	main.roc <- if(is.null(main)) {
        switch(class(x)[1], "AROC.kernel" = "AROC Kernel-based", "AROC.bnp" = "AROC Bayesian nonparametric", "AROC.sp" = "AROC semiparametric")
    } else {
        main
    }
	main.auc <- "AAUC"

    type.lines <- if(class(x)[2] %in% c("AROC.kernel", "AROC.sp")) {
        "s"
    } else {
        "l"
    }
    
    plot(x$p, x$ROC[,1], xlab = "FPF", ylab = "TPF", xlim = c(0,1), ylim = c(0,1), main = main.roc, type = type.lines, cex.lab = 1.3, cex.axis = 1.3, ...)
    if(ncol(x$ROC) == 3) {
        lines(x$p, x$ROC[,2], lty = 2, type = type.lines)
        lines(x$p, x$ROC[,3], lty = 2, type = type.lines)
    }
	abline(0,1, col = "grey")
	if(length(x$AUC) == 3) {
		legend.text <- paste0(main.auc, ": ", paste(round(x$AUC[1], 3), " (", round(x$AUC[2], 3),"",", ", round(x$AUC[3], 3),")", sep = ""))
	} else {
		legend.text <- paste0(main.auc, ": ", round(x$AUC[1], 3))
	}
	legend(0.4, 0.2, legend.text, bty = "n", cex = if(!is.null(dots$cex)) dots$cex else 1.3)

}
