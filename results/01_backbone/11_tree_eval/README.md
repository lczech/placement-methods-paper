Overview
-------------------------

Here, we provide code for the evaluations of the automatic reference trees.
We used EPA-ng to place the original SILVA sequences from which the trees were created on those trees.
The trees were created using the consensus sequences of sequences at the tips of the underlying taxonomy.
Thus, if the trees are superduper nice, we'd expect all sequences 
to be placed on the exact same branch into whose consensus sequence they went.

Placement
-------------------------

The sequences of Silva are prepared using the program `code/art/silva_tree_eval_queries.cpp`,
in order to get into a format and naming convention that is easy to use for our
evaluation program later.

We used the EPA-ng binary dump feature for speedup,
thus, first the `binary_dump.sh` script has to be run.
Then, use the `job_simple.sh` script for small datasets like the `Archaea` to place the Silva sequences.
For the larger datasets, the `create_jobs.sh` script uses `job_template.sh` to
create smaller jobs, so that the run time is reduced.
These produce placement files to be used for our evaluation.

Evaluation
-------------------------

The evaluation is done by using the program
`code/art/silva_tree_eval.cpp`

This creates log scripts in subdirectories for each tree that contain a lot of detailled data.
Also it creates our needed data in the `viz` dir:

`viz/lists`: lists of distances (in num of edges and branch length units) for all placements.
`vis/tables`: summarizes the lists into histograms. not longer needed, as this is done in python now.

These files are then used to make the plots.

Visualization
-------------------------

Now, it's time to produce the graphs of the paper, showing
the placement accuracy, that is, measuring the distance of the Silva sequences
from their expected placement branches.

For this, we used the `make_hists.py` script to turn the lists into histograms for the paper.
The figures end up in a subdirectory `figures`.

Misc
-------------------------

Earlier tests used the `eval_all.sh` (and the `log`) to run the eval individually for the directories.
Now, this is all done within the `silva_tree_eval` itself.


Commands
-------------------------

Merge `jplace` files into one:

`./merge_jplace_files job_*/epa_result.jplace epa_result.jplace > merge_jplace_files.log`

Make evaluation statistics (old):

`./silva_tree_eval &> eval_full.log`

Make evaluation statistics (old):

`./silva_tree_eval epa_result.jplace Archaea > silva_tree_eval.log`

Make placement histograms:

`./placement_histograms epa_result.jplace /path/to/data/`

