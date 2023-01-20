## Resubmission
This is a resubmission. In this version I have:

* Removed the \dontrun{} calls from all examples. When checking the 
package at https://win-builder.r-project.org/, two examples, in  
R/find_bayes_estimate.R and R/multilink-package.R, took longer than 5 sec to
execute without \dontrun or \donttest calls. I have thus wrapped portions of 
these examples with \donttest{} calls. When rechecking the 
package at https://win-builder.r-project.org/ with these changes all examples
executed in < 5 sec. If any examples take longer than 5 sec to execute on this 
submission I am happy to shorten the execution times and resubmit.

* Added verbose arguments to R/create_comparison_data.R,
R/find_bayes_estimate.R, and R/gibbs_sampler.R, that wrap all print statements
in these functions. When verbose is set to FALSE, all print statements are 
suppressed for the given function. I have set the default for the verbose 
arguments to TRUE, as these functions can take some time to run for moderate 
sized data sets, and thus progress messages are important for users to know when 
to expect the functions to finish running.

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
