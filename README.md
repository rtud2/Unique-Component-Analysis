# Unique Component Analysis (UCA)

A technique to carry out contrastive Principal Component Analysis with multiple backgrounds based on Lagrange Multipliers.
UCA can be used with multiple background data to remove known variation from multiple sources and implemented to handle high-dimensional data. 

##  Download
You may compile directly from source via command line or shell, using the UCA folder with `R CMD build UCA`, or download the pre-compiled files `uca_0.13.zip` for windows or `uca_0.13.tar.gz` for linux machines. To install the precompiled versions,

```
install.packages("uca_0.13.zip", repos=NULL, type = "source") #replace with uca_0.13.tar.gz if on linux
library('uca')
```

The package depends on `Rfast`, `Rcpp`, and `RSpectra`.

## functions
Helper functions useful for analyzing image KDEF faces examples.

## examples
The example folder has the analysis corresponding to the paper. 


