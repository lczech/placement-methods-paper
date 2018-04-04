Overview
-------------------------

This directory contains the code and scripts that we used for the papers.
That is, using this, one should be able to reproduce our exact results
as presented in the papers.
We modified our original code by removing hard coded file paths and
similar details that would have made reproduction of results harder.
By this, we tried to make the code as portable as feasible.
However, it is still possible that some passages have to be adjusted
for other environments in order to work properly.

`art`
-------------------------

Our prototype implementation of the Automatic Reference Tree method,
as well as the code that we used for the testing and evaluation.

`data`
-------------------------

Data preprocessing programs that are not specific to one method,
but were used to get the sequence data into usable formats in the first place.

`tests`
-------------------------

Some old code and scripts that we used during the development process.
These are just provided as-is, and probably won't compile without some changes.
We keep them here just in case we later need to revisit some of the
early details of our methods.

`visualization`
-------------------------

Prototypes of the visualization methods, e.g, Edge Dispersion and Edge Correlation,
as well as related programs.
