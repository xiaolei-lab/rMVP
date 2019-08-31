# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


logging.initialize <- function(module) {
    file <- NULL
    if (options("rMVP.OutputLog2File") == TRUE) {
        now <- Sys.time()
        file <- paste(module, format(now, "%Y%m%d_%H%M%S"), "log", sep = ".")
    }
    
    assign("logging.file", file, envir = package.env)
}


logging.log <- function(..., file = NULL, sep = " ", fill = FALSE, labels = NULL, verbose = TRUE) {
    if (verbose) {
        cat(..., sep = sep, fill = fill, labels = labels)
    }
    
    if (is.null(file)) {
        try(file <- get("logging.file", envir = package.env), silent = TRUE)
    }
    
    if (!is.null(file)) {
        cat(..., file = file, sep = sep, fill = fill, labels = labels, append = TRUE)
    }
}