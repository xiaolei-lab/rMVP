.onLoad <- function(libname, pkgname) {
    op <- options()
    op.rMVP <- list(
        rMVP.OutputLog2File = TRUE
    )
    toset <- !(names(op.rMVP) %in% names(op))
    if (any(toset)) { 
        options(op.dplyr[toset])
    }
    
    return(invisible())
}

