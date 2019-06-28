.onLoad <- function(libname, pkgname) {
    # Limit number of threads in veclib (MacOS MRO)
    if (Sys.info()["sysname"] == "Darwin") {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    
    op <- options()
    op.rMVP <- list(
        rMVP.OutputLog2File = TRUE
    )
    toset <- !(names(op.rMVP) %in% names(op))
    if (any(toset)) { 
        options(op.rMVP[toset])
    }
    
    
    return(invisible())
}
