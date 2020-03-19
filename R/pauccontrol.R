pauccontrol <-
function(compute = FALSE, focus = c("FPF", "TPF"), value = 1)
  list(compute = compute, focus = match.arg(focus), value = value)
