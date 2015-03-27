getSymbolicNamesList <- function(strings) {
  lapply(strings, function(v) deparse(as.name(v), backtick=TRUE))
}
