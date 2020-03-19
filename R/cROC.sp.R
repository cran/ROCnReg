cROC.sp <-
function(formula.healthy, formula.diseased, group, tag.healthy, data,  newdata, est.cdf = c("normal", "empirical"), pauc = pauccontrol(), p = seq(0,1,l = 101), B = 1000) {
    
    pauc <- do.call("pauccontrol", pauc)
    
    compute.ROC <- function(formula.healthy, formula.diseased, group, tag.healthy, data, newdata, est.cdf, pauc, p = seq(0,1,l = 101)) {
        data.h <- data[data[,group] == tag.healthy,]
        data.d <- data[data[,group] != tag.healthy,]
        
        n0 <- nrow(data.h)
        n1 <- nrow(data.d)
        
        np <- length(p)
        
        marker <- all.vars(formula.healthy)[1]
        
        # Fit the model in the healthy population
        fit0p <- lm(formula = formula.healthy, data = data.h)
        sigma0p <- summary(fit0p)$sigma
        
        
        # Fit the model in the diseased population
        fit1p <- lm(formula = formula.diseased, data = data.d)
        sigma1p <- summary(fit1p)$sigma

        # Coefficients ROC curve
		coeffs <- c(names(coefficients(fit0p)), names(coefficients(fit1p))[is.na(match(names(coefficients(fit1p)), names(coefficients(fit0p))))])
		beta.h <- beta.d <- rep(0, length(coeffs))
		names(beta.h) <- names(beta.d) <- coeffs
		
		beta.h[match(names(coefficients(fit0p)), coeffs)] <- coefficients(fit0p)
		beta.d[match(names(coefficients(fit1p)), coeffs)] <- coefficients(fit1p)

		beta.ROC <- c((beta.h - beta.d)/sigma1p, "b" = sigma0p/sigma1p)
        
        # Three terms of the cROC
        a <- (predict(fit0p, newdata = newdata) - predict(fit1p, newdata = newdata))/sigma1p
        b <- sigma0p/sigma1p
        
        if(est.cdf == "normal") {
            cROC <- 1 - apply(a + outer(rep(b, nrow(newdata)), qnorm(1-p), "*"), c(1,2), pnorm)
            #cAUC <- apply(cROC, 1, simpson, set.p = p)
            # Binormal model
            cAUC <- 1 - pnorm(a/sqrt(1+b^2))
			if(pauc$compute){
				if(pauc$focus == "FPF"){
					#pu <- seq(0, pauc$value, len = np)
					#cROC_pauc <- 1 - apply(a + outer(rep(b, nrow(newdata)), qnorm(1-pu), "*"), c(1,2), pnorm)
					#cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
                    cpAUC <- pbivnorm(-a/sqrt(1+b^2), qnorm(pauc$value), -b/sqrt(1+b^2))
				} else{
					#a <- (predict(fit1p, newdata = newdata) - predict(fit0p, newdata = newdata))/sigma0p
					#b <- sigma1p/sigma0p
					#pu <- seq(pauc$value, 1, len = np)
					#cROC_pauc <- apply(a + outer(rep(b, nrow(newdata)), qnorm(1-pu), "*"), c(1,2), pnorm)
					#cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
                    cpAUC <- pbivnorm(-a/sqrt(1+b^2), qnorm(1-pauc$value), -1/sqrt(1+b^2))
				}
			}

        } else {
            #res0p <- fit0p$residuals/sigma0p
            #F0res <- ecdf(res0p)
            
            #res1p <- fit1p$residuals/sigma1p
            #F1res <- ecdf(res1p)
            
            #cdf0 <- apply(outer(res0p, res0p, "<="), 2, mean)
            #cdf0_inv <- apply(outer(cdf0, 1-p, ">="), 2, function(x, z) {
            #        res <- min(z[x]) #min(c(z[x], max(z)))
            #        res
            #    }, z = res0p)
            #cdf0_inv <- replace(cdf0_inv, is.infinite(cdf0_inv), max(res0p))
            
            res0p <- fit0p$residuals/sigma0p
            res1p <- fit1p$residuals/sigma1p
            F1res <- ecdf(res1p)
            cdf0_inv <- quantile(res0p, 1-p, type = 1)
            
            # cROC: rows covariates / columns FPF
            cROC <- 1 - apply(a + outer(rep(b, nrow(newdata)), cdf0_inv, "*"), c(1,2), F1res)
            cAUC <- apply(cROC, 1, simpson, set.p = p)
			if(pauc$compute){
				if(pauc$focus =="FPF"){
					pu <- seq(0, pauc$value, len = np)
					cdf0_inv_pauc <- quantile(res0p, 1-pu, type = 1)
					cROC_pauc <- 1 - apply(a + outer(rep(b, nrow(newdata)), cdf0_inv_pauc, "*"), c(1,2), F1res)
					cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
				} else{
					a <- (predict(fit1p, newdata = newdata) - predict(fit0p, newdata = newdata))/sigma0p
					b <- sigma1p/sigma0p
					pu <- seq(pauc$value, 1, len = np)
					cdf1_inv_pauc <- quantile(res1p, 1-pu, type = 1)
					F0res <- ecdf(res0p)
					cROC_pauc <- apply(a + outer(rep(b, nrow(newdata)), cdf1_inv_pauc, "*"), c(1,2), F0res)
					cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
				}
			}
        }
        res <- list()
        res$p <- p
        res$ROC <- cROC
        res$AUC <- cAUC
        if(pauc$compute){   
            res$pAUC <- cpAUC
        }
        res$data.h <- data.h
        res$data.d <- data.d
        res$fit <- list(h = fit0p, d = fit1p)
        res$coeff <- list(h = coefficients(fit0p), d = coefficients(fit1p), ROC = beta.ROC)
        res
    }
    
    est.cdf <- match.arg(est.cdf)
    np <- length(p)
    
    # Check if the seq of FPF is correct
    if(est.cdf == "empirical") {
        if(min(p) != 0 | max(p) != 1 | np%%2 == 0) {
            warning("The set of FPFs at which the covariate-specific ROC curve is estimated is not correct. The set is used for calculating the AUC using Simpson's rule. As such, it should (a) start in 0; (b) end in 1; and (c) have an odd length.")
        }
    }

    if(inherits(formula.healthy, "character")) {
        formula.healthy <- as.formula(formula.healthy)
    }
    if(inherits(formula.diseased, "character")) {
        formula.diseased <- as.formula(formula.diseased)
    }
    # Marker variable
    tf <- terms.formula(formula.healthy, specials = c("f"))
    if (attr(tf, "response") > 0) {
        marker.h <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The 'formula.healthy' should include the response variable (left hand side)")
    }
    
    tf <- terms.formula(formula.diseased, specials = c("f"))
    if (attr(tf, "response") > 0) {
        marker.d <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The 'formula.healthy' should include the response variable (left hand side)")
    }
    
    marker <- marker.h
    
    # Variables in the model
    # Variables in the model
    names.cov.h <- all.vars(formula.healthy)[-1]
    names.cov.d <- all.vars(formula.diseased)[-1]
    names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])
    
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
    	stop("Newdata must be a data frame")
    if(sum(is.na(match(c(marker, names.cov, group), names(data)))))
    	stop("Not all needed variables are supplied in data")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
    	stop("Not all needed variables are supplied in newdata")
    if(length(unique(data[,group])) != 2)
    	stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep = "")
    
    # Data, removing missing values
    data.new <- data[,c(marker,group,names.cov)]
    omit.h <- apply(data.new[data.new[,group] == tag.healthy, c(marker, group, names.cov)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.healthy, c(marker, group, names.cov)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.healthy,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.healthy,,drop = FALSE][!omit.d,,drop = FALSE])
    
    # New data (for predictions)
    if(missing(newdata)) {
        newdata <- cROCData(data.new, names.cov, group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop=FALSE])
    }
    
    res.fit <- compute.ROC(formula.healthy = formula.healthy, formula.diseased = formula.diseased, group = group, tag.healthy = tag.healthy, data = data.new, newdata  = newdata, est.cdf = est.cdf, pauc = pauc, p = p)
    crocp <- res.fit$ROC
    caucp <- res.fit$AUC
    if(pauc$compute){
        cpaucp <- res.fit$pAUC
    }
    coeff0p <- res.fit$coeff$h
    coeff1p <- res.fit$coeff$d
    coeffROCp <- res.fit$coeff$ROC

    if(B > 0) {
        # Confidence intervals
        crocpb <- array(0, dim = c(nrow(newdata), length(p), B))
        caucb <- matrix(0, nrow = nrow(newdata), ncol = B)
        if(pauc$compute){
            cpaucb <- matrix(0, nrow = nrow(newdata), ncol = B)
        }
        ccoeff0b <- matrix(0, nrow = length(coeff0p), ncol = B)
        ccoeff1b <- matrix(0, nrow = length(coeff1p), ncol = B)
        ccoeffROCb <- matrix(0, nrow = length(coeffROCp), ncol = B)
        for(l in 1:B) {
            data.boot.d <- res.fit$data.d
            data.boot.h <- res.fit$data.h
            
            res.h.b <- sample(res.fit$fit$h$residuals, replace = TRUE)
            res.d.b <- sample(res.fit$fit$d$residuals, replace = TRUE)
            
            data.boot.h[,marker] <- res.fit$fit$h$fitted + res.h.b
            data.boot.d[,marker] <- res.fit$fit$d$fitted + res.d.b
            
            data.boot <- rbind(data.boot.d, data.boot.h)
            
            res.boot <- compute.ROC(formula.healthy = formula.healthy, formula.diseased = formula.diseased, group = group, tag.healthy = tag.healthy, data = data.boot, newdata = newdata, est.cdf = est.cdf, pauc = pauc, p = p)
            crocpb[,,l] <- res.boot$ROC
            caucb[,l] <- res.boot$AUC
            if(pauc$compute){
                cpaucb[,l] <- res.boot$pAUC
            }
            ccoeff0b[,l] <- 	res.boot$coeff$h
    		ccoeff1b[,l] <- res.boot$coeff$d
    		ccoeffROCb[,l] <- res.boot$coeff$ROC
        }
    }
    
    columns <-switch(as.character(B > 0),"TRUE" = 1:3, "FALSE" = 1)
    col.names <-c("AUC","AUCql", "AUCqh")[columns]
    col.names_coeff <-c("est","ql", "qh")[columns]
    if(pauc$compute){
        col.names_pauc <-c("pAUC","pAUCql", "pAUCqh")[columns]
    }
    cROCres <- list()
    cROCres$est <- crocp
    cAUCres <- matrix(0, ncol = length(columns), nrow = nrow(newdata), dimnames = list(1:nrow(newdata), col.names))
    cAUCres[,1] <- caucp
    if(pauc$compute){
        cpAUCres <- matrix(0, ncol = length(columns), nrow = nrow(newdata), dimnames = list(1:nrow(newdata), col.names))
        cpAUCres[,1] <- cpaucp
    }
    ccoeff0res <- matrix(0, ncol = length(columns), nrow = length(coeff0p), dimnames = list(names(coeff0p), col.names_coeff))
    ccoeff0res[,1] <- coeff0p

    ccoeff1res <- matrix(0, ncol = length(columns), nrow = length(coeff1p), dimnames = list(names(coeff1p), col.names_coeff))
    ccoeff1res[,1] <- coeff1p

    ccoeffROCres <- matrix(0, ncol = length(columns), nrow = length(coeffROCp), dimnames = list(names(coeffROCp), col.names_coeff))
    ccoeffROCres[,1] <- coeffROCp

    if(B > 0) {
        cROCres$ql <- apply(crocpb, c(1,2), quantile, prob = 0.025)
        cROCres$qh <- apply(crocpb, c(1,2), quantile, prob = 0.975)
        
        cAUCres[,2] <- apply(caucb, 1, quantile, prob = 0.025)
        cAUCres[,3] <- apply(caucb, 1, quantile, prob = 0.975)
        
        if(pauc$compute){
            cpAUCres[,2] <- apply(cpaucb, 1, quantile, prob = 0.025)
            cpAUCres[,3] <- apply(cpaucb, 1, quantile, prob = 0.975)
        }

        ccoeff0res[,2] <- apply(ccoeff0b, 1, quantile, prob = 0.025)
        ccoeff0res[,3] <- apply(ccoeff0b, 1, quantile, prob = 0.975)

        ccoeff1res[,2] <- apply(ccoeff1b, 1, quantile, prob = 0.025)
        ccoeff1res[,3] <- apply(ccoeff1b, 1, quantile, prob = 0.975)

        ccoeffROCres[,2] <- apply(ccoeffROCb, 1, quantile, prob = 0.025)
        ccoeffROCres[,3] <- apply(ccoeffROCb, 1, quantile, prob = 0.975)
    }
    colnames(cAUCres) <- col.names
    if(pauc$compute){
        colnames(cpAUCres) <- col.names_pauc
    }
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$data <- data
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$marker <- marker
    res$group <- group
    res$tag.healthy <- tag.healthy
    res$formula <- list(h = formula.healthy, d = formula.diseased)
    res$est.cdf <- est.cdf
    res$ci.fit <- ifelse(B, TRUE, FALSE)
    res$p <- p
    res$ROC <- cROCres
    res$AUC <- cAUCres
    if(pauc$compute){
        res$pAUC <- if(pauc$focus == "FPF") {
        	cpAUCres/pauc$value
        } else {
        	cpAUCres/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    res$fit <- res.fit$fit
    res$coeff <- list(h = ccoeff0res, d = ccoeff1res, ROC = ccoeffROCres)
    class(res) <- c("cROC","cROC.sp")
    res
}