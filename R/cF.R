cF <-
function(x, iter, Xpred, icov, P, Beta, Sigma) {
	res <- sum(P[iter,]*pnorm(x, mean = Xpred[icov,]%*%t(Beta[iter,,]), sd = Sigma[iter,]))
	res
}
