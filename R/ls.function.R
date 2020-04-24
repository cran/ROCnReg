ols.function <- function(X, y, vcov = FALSE) {
	res <- list()
	if(vcov) {
		res$vcov <- solve(crossprod(X)) 
		res$coeff <- res$vcov%*%crossprod(X,y)
	} else {
		res$coeff <- try(solve(crossprod(X), crossprod(X,y)), silent = TRUE)	
	}
	res	
}
