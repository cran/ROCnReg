qFk <-
function(p, y, h){
	toInvert <- function(x, p, y, h) {
	  return(mean(pnorm((x-y)/h))-p)
	}
	res <- uniroot(toInvert, interval = c(min(y) - 10^25, max(y) + 10^25), p, y, h)$root
	return(res) 
}
