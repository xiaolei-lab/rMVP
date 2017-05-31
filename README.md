# MVP
## A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study
### Installation
**MVP** is not available on CRAN, but can be installed using the R package **devtools**. For installing **MVP**, some packages needed to be installed beforehand. **MVP** can be installed with the following R code:
```r
source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("xiaoleiLiubio/MVP")
```
After installed successfully, **MVP** can be loaded with the following R code:
```r
library(MVP)
```
