qF <-
function(q, iter, Xpred, icov, P, Beta, Sigma) {
	toInvert <-  function(x, q, iter, Xpred, icov, P, Beta, Sigma) {
		return(sum(P[iter,]*pnorm(x, mean = Xpred[icov,]%*%t(Beta[iter,,]), sd = Sigma[iter,])) - q)
	}
	res <-  uniroot(toInvert, interval = c(-10^15,10^15), q, iter, Xpred, icov, P, Beta, Sigma)$root
	return(res)
}
