package.env <- NULL

.onLoad <- function(libname, pkgname) {
    # Limit number of threads in veclib (MacOS MRO)
    if (Sys.info()["sysname"] == "Darwin") {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    } 
    
    # set option
    op <- options()
    op.rMVP <- list(
        rMVP.OutputLog2File = TRUE
    )
    toset <- !(names(op.rMVP) %in% names(op))
    if (any(toset)) { 
        options(op.rMVP[toset])
    }
    
    # package level environment
    package.env <<- new.env()
    
    return(invisible())
}

.onAttach <- function(...){
    packageStartupMessage("Full description, Bug report, Suggestion and the latest version:")
    packageStartupMessage("https://github.com/xiaolei-lab/rMVP")
}
