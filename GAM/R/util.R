setdiff.data.frame <-
    function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]

#' Load data only if its absent from global environment
#' @param name Name of a dataset
#' @param ... Additional arguments for data()
lazyData <- function(name, ...) {
    if (!name %in% ls(envir=.GlobalEnv)) {
        print(paste0("No ", name, ", loading"))
        do.call(data, list(name, ...))
    }
}
