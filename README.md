# Unique Component Analysis (UCA)

A technique to carry out contrastive Principal Component Analysis with multiple backgrounds based on Lagrange Multipliers.
UCA can be used with multiple background data to remove known variation from multiple sources and implemented to handle high-dimensional data. 

##  Download
You may compile directly from source via command line or shell, using the UCA folder with `R CMD build UCA`. Note, The package depends on `Rcpp`, `RcppArmadillo`, and `RSpectra`, so you will need those installed before you can compile.

 You may also choose to download the pre-compiled files `uca_0.14.zip` for windows or `uca_0.14.tar.gz` for linux machines. To install the precompiled versions,

```
install.packages("uca_0.13.zip", repos=NULL, type = "source") #replace with uca_0.13.tar.gz if on linux
library('uca')
```

## functions
Helper functions useful for analyzing image KDEF faces examples.

## examples
The example folder has the analysis corresponding to the paper. The mouse folder contains the relevant UCA code for our results in our paper.

The mouse protein expression dataset can be downloaded at [https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression](https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression).

The KDEF data used in the faces analysis can be downloaded at [https://www.kdef.se/](https://www.kdef.se/)

The timing folder contains a simulation that demonstrates the speed increase using the Product SVD method, rather than the covariance method, to run `UCA`.

## Changelogs
2/1/21 - remove dependency on Rfast
