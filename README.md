# MVP
## A Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool For Genome-Wide Association Study
### Installation
**MVP** is not available on CRAN, but can be installed using the R package **"devtools"**. There are two packages should be installed beforehand, **"snpStats"** and **"rfunctions"**. **MVP** can be installed with the following R code:
```r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_github("Bioconductor-mirror/snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("xiaoleiLiubio/MVP")
```
After installed successfully, **MVP** can be loaded with the following R code:
```r
library(MVP)
```
