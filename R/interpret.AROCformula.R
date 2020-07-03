interpret.AROCformula <-
function(formula, data) {
	env <- environment(formula)
    encl <- loadNamespace("ROCnReg")

    if(inherits(formula, "character"))        
        formula <- as.formula(formula)

    tf <- terms.formula(formula, specials = c("f", "ns", "bs"))
    if(!is.null(attr(tf,"specials")$ns) | !is.null(attr(tf,"specials")$bs)) {
        stop("'ns' (natural splines) or 'bs' (B-splines) are not allowed in the formula. Please use 'f' instead.")
    } 

    if (attr(tf, "response") > 0) {
        marker <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The formula should include the response variable (left hand side)")
    }

    terms <- attr(tf, "term.labels")
    nt <- length(terms)  
    
    ifun <- sort(attr(tf,"specials")$f) - 1
    nfun <- length(ifun)
    smooth <- terms[ifun]

    if(nfun > 0) {
        fixed <- terms[-ifun]
        ilin <- (1:nt)[-ifun]
        nlin <- length(fixed)
    } else {
        fixed <- terms
        ilin <-  1:nt
        nlin <- length(fixed)
    }

    II <- list()
    h  <- list()
    K_total <- list()
    partial <- vector()
    k <- 0
    vars.formula <- NULL
    if(nt) {        
        if(nfun > 0) {
            for(i in ifun) {
                k <- k + 1                   
                st <- eval(parse(text = paste("AROC.", terms[i], sep = "")))
                II[[k]] <- st$cov
                h[[k]] <- -1
                K_total[[k]] <- st$K
                partial[k] <- terms[i]
                vars.formula <- c(vars.formula, st$term)           
            }
        }
        # Parametric (linear and categorical: all in one)
        if(nlin > 0) {
            k <- k + 1
            full_term <- paste(terms[ilin], collapse = "+", sep = "")
            II[[k]]<- c(-1, full_term)
            h[[k]] <- 0 # parametric
            K_total[[k]] <- 0
            partial[k] <- full_term
            vars.formula <- c(vars.formula, all.vars(formula(paste("~", full_term))))

        }
        #names.cov <- all.vars(formula)[-1]
        names.cov <- vars.formula
        data.cov <- data[, names(data) %in% names.cov, drop = FALSE]
        #numeric.var <- names.cov[!apply(data.cov, 2, is.factor)]
        numeric.var <- names.cov[!sapply(names.cov, function(x, data) is.factor(data[,x]), data = data.cov)]
        if(length(numeric.var) != 0) {
            cov.std <- matrix(ncol = length(numeric.var), nrow = 2, dimnames = list(c("Mean", "sd"), numeric.var))
            data.cov.std <- data.cov
            for(i in 1:length(numeric.var)) {
                cov.std[1,i] <- mean(data.cov[,numeric.var[i]], na.rm = TRUE)
                cov.std[2,i] <- sd(data.cov[,numeric.var[i]], na.rm = TRUE)
                data.cov.std[,numeric.var[i]] <- (data.cov[,numeric.var[i]] - cov.std[1,i])/cov.std[2,i]
            }
        } else {
            cov.std <- NULL
            data.cov.std <- data.cov
        }
    } else { # Only the intercept
        names.cov <- NULL
        data.cov <- NULL
        cov.std <- NULL
        data.cov.std <- NULL
    }

    II <- if(length(II)) {
        matrix(unlist(II), nrow = 2)
    } else {
        matrix(0, nrow = 2)
    }      
    #res <- list(marker = marker, II = II, h = unlist(h), K = unlist(K), npartial = k, partial = partial, data.cov = data.cov, cov.std = cov.std, data.cov.std = data.cov.std)
    res <- list(marker = marker, II = II, h = unlist(h), K = K_total, npartial = k, partial = partial, data.cov = data.cov, cov.std = cov.std, data.cov.std = data.cov.std)
    res
}
