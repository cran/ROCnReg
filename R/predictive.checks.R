predictive.checks <-
function(object, statistics = c("min","max","kurtosis","skewness"), ndensity = 512, devnew = TRUE) {
	UseMethod("predictive.checks")
}
