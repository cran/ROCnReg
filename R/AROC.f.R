AROC.f <-
function(..., by = NULL, K = 0) {
	vars <- as.list(substitute(list(...)))[-1] # List
	args <- match.call()

	d <- length(vars)
	if(d == 0 | d > 1) {
		stop("Incorrect number of covariates")
	}

	term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
	if (term[1] == ".") {
		stop("f(.) not yet supported.")
	}

	for (i in 1:d){
		term[i] <- attr(terms(reformulate(term[i])), "term.labels")
	}	

	#if(is.null(args$x1) & is.null(args$x2))
	#	stop("x1 or x2 must be indicated")
	#if(!is.null(args$x1) & is.null(args$by)) { # Smooth effect			   
	#	cov = c("-1", deparse(args$x1, backtick = TRUE, width.cutoff = 500))
	#} else if (!is.null(args$x1) & !is.null(args$by)) {	  
	#	cov = c(deparse(args$by, backtick = TRUE, width.cutoff = 500), deparse(args$x1, backtick = TRUE, width.cutoff = 500))
	#} else {
	#	stop("Invalid expression")
	#}

	if(d == 1 & is.null(args$by)) { # 1D Smooth effect			   
		cov <- c("-1", term[1])
	} else if (d == 1 & !is.null(args$by)) { # Factor by curve 	  
		cov <- c(deparse(args$by, backtick = TRUE, width.cutoff = 500), term[1])
	} else {
		stop("Invalid expression")
	}

	res <- list(term = term, cov = cov, K = K)
	res
}
