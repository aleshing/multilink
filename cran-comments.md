## Resubmission Notes
Thank you for reviewing my package submission. The review of my previous 
submission was as follows:

"\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary.
Please replace \dontrun with \donttest.

Please unwrap the examples if they are executable in < 5 sec, or replace
dontrun{} with \donttest{}.

You write information messages to the console that cannot be easily
suppressed.
It is more R like to generate objects that can be used to extract the
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning() or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console. (except for print, summary, interactive functions)
-> R/create_comparison_data.R; R/find_bayes_estimate.R; R/gibbs_sampler.R

Please fix and resubmit."

I have addressed these issues as follows:

* I have removed the \dontrun{} calls from all examples. When checking the 
package at https://win-builder.r-project.org/, two examples, in  
R/find_bayes_estimate.R and R/multilink-package.R, took longer than 5 sec to
execute without \dontrun or \donttest calls. I have thus wrapped portions of 
these examples with \donttest{} calls. When rechecking the 
package at https://win-builder.r-project.org/ with these changes all examples
executed in < 5 sec. If any examples take longer than 5 sec to execute on this 
submission I am happy to shorten the execution times and resubmit.

* I have added verbose arguments to R/create_comparison_data.R,
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
