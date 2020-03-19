predict.design.matrix.bnp <-
function(object, newdata, ...) {
    if(object$iformula$npartial == 0) { # Only the intercept
        Xp <- matrix(1, ncol = 1, nrow = nrow(newdata))
    } else {
        Xp <- NULL    
        # Organize the newdataframe as it was in the original data
        cov.names <- names(object$iformula$data.cov)
        newdata <- newdata[, cov.names, drop = FALSE]

        # Standardised the continuous covariates
        newdata.std <- newdata
        cov.names.std <- colnames(object$iformula$cov.std)
        if(!is.null(cov.names.std)) {
            for(i in 1:length(cov.names.std)) {
                aux <- object$iformula$cov.std[,cov.names.std[i]]
                newdata.std[, cov.names.std[i]] <- (newdata[,cov.names.std[i]] - aux[1])/aux[2]
            }
        }
        for(i in 1:object$iformula$npartial) {
            if(any(object$iformula$II[,i] == -1)) {
                if(object$iformula$h[i] == 0 | object$iformula$h[i] == 1) { # Linear and factor
                    if(object$standardise) {
                        mfp <- model.frame(object$terms[[i]], newdata.std, xlev = attr(object$terms[[i]], "xlev"))
                    } else {
                        mfp <- model.frame(object$terms[[i]], newdata, xlev = attr(object$terms[[i]], "xlev"))
                    }
                    Xp_aux <- model.matrix(object$terms[[i]], data = mfp, contrasts.arg = attr(object$terms[[i]], "contrast"))[,-1, drop = FALSE]
                    Xp <- cbind(Xp, Xp_aux)
                } else {
                    Bs <- suppressWarnings(predict.bbase.os(object$terms[[i]], newdata[,object$iformula$II[2,i]]))
                    Xp <- cbind(Xp, Bs)
                }
            } else { # Factor by curve
                Bs <- predict.bbase.interaction.factor.by.curve.os(object$terms[[i]], newdata[,object$iformula$II[2,i]], newdata[,object$iformula$II[1,i]])
                Xp <- cbind(Xp, Bs)
            }
        }
        # Add the intercept
        Xp <- cbind(1, Xp)
    }
    res <- list()
    res$X <- Xp
    res
}
