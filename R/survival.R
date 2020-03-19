survival <-
function(y, obs) {
    n <- length(y)
    estimate <- (sum(y > obs) + 0.5*sum(y == obs))/n
}
