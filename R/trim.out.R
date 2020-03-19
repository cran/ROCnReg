trim.out <-  function (x, trim = 0.1, na.rm = FALSE) {
    if (na.rm) 
        x <- x[!is.na(x)]
    if (!is.numeric(trim) || length(trim) != 1L) 
        stop("'trim' must be numeric of length one")
    n <- length(x)
    if (trim > 0 && n) {
        if (is.complex(x)) 
            stop("trim is not defined for complex data")
        if (anyNA(x)) 
            return(NA_real_)
        if (trim >= 0.5 && trim < 1) 
            return(NA_real_)
        if (trim < 1) 
            lo <- floor(n * trim) + 1
        else {
            lo <- trim + 1
            if (trim >= (n/2)) 
                return(NA_real_)
        }
        hi <- n + 1 - lo
        res <- sort.int(x, index.return = TRUE)
        trimi <- res[["ix"]][c(1:(lo - 1), (hi + 1):length(x))]
        x <- res[["x"]][lo:hi][order(res[["ix"]][lo:hi])]
        attr(x, "trim") <- trimi
    }
    return(x)
}