## Resubmission
This is a resubmission. In this version I have fixed the clang-UBSAN related
additional issue that arose during CRAN package checks. 

As described in an email to the CRAN team, on the MacOS architectures, the CHECK 
returns a NOTE because the libs subdirectory is then above the 1MB threshold. 
However, it seems that this NOTE only appears under MacOS but not under Windows 
or Linux. My understanding is that the libs subdirectory is large due to the 
use of Rcpp. Some functions of the multilink package have been written in C++ 
using Rcpp to perform MCMC. Without the speed up gained from those C++ 
functions, this package would become impractical. 

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
