get_vars_formula <- function(formula) {
	env <- environment(formula) 
    if(inherits(formula, "character"))        
        formula <- as.formula(formula)

    tf <- terms.formula(formula, specials = c("f", "ns", "bs"))
    #if(!is.null(attr(tf,"specials")$ns) | !is.null(attr(tf,"specials")$bs)) {
    #    stop("'ns' (natural splines) or 'bs' (B-splines) are not allowed in the formula. Please use 'f' instead.")
    #}

    if (attr(tf, "response") > 0) {
        marker <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The formula should include the response variable (left hand side)")
    }

    # All terms
    terms <- attr(tf, "term.labels")
    nt <- length(terms)  
    
    # Smooth terms
    ifun <- sort(attr(tf,"specials")$f) - 1
    nfun <- length(ifun)

    if(nfun > 0) {
        fixed <- terms[-ifun]
        ilin <- (1:nt)[-ifun]
        nlin <- length(fixed)
    } else {
        fixed <- terms
        ilin <-  1:nt
        nlin <- length(fixed)
    }
	vars.formula <- NULL
    if(nt) {        
        if(nfun > 0) {
            for(i in ifun) {
                vars.formula <- c(vars.formula, eval(parse(text = paste("AROC.", terms[i], sep="")))$term)   
            }
		}
		if(nlin > 0) {
            full_term <- paste(terms[ilin], collapse = "+", sep = "")
            vars.formula <- c(vars.formula, all.vars(formula(paste("~", full_term))))
        }
    }
    vars.formula
}