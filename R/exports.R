#' @rawNamespace exportPattern("^MVP")
NULL

#' @useDynLib rMVP, .registration=TRUE
NULL
#> NULL

#' @import methods
NULL

#' @import stats
NULL

#' @import MASS
NULL

#' @import graphics
NULL

#' @import grDevices
NULL

#' @import bigmemory
NULL

#' @importFrom parallel detectCores
NULL

#' @importFrom utils read.delim read.table write.table write.csv memory.limit packageVersion
NULL

#' @importFrom bigmemory is.big.matrix
#' @export
bigmemory::is.big.matrix

#' @importFrom bigmemory as.big.matrix
#' @export
bigmemory::as.big.matrix

#' @importFrom bigmemory attach.big.matrix
#' @export
bigmemory::attach.big.matrix

#' @importFrom bigmemory deepcopy
#' @export
bigmemory::deepcopy

#' @importFrom methods new
#' @export
methods::new
