Overview
-------------------------

This directory contains the script used to generate per-sample placement files,
using the chunk and batch files of the previous steps.
It uses the prototype program `jplace_unchunkify`.

As we used an experiment version of EPA-ng, we dicovered a bug in the output jplace files.
Thus, the process script is a bit complicated in order to remove the effects of this bug.
The placement results themselves were not affected by this,
except for a few outlier sequences (which caused the bug to appear in the first place,
so excluding them did not change anything - their placement was bogus anyway).

Once processed, the directory for each dataset contains the following stuff:

 * `check/`		intermediate files for checking whether the created sample jplace files contain all needed entries (that is, all from the map, except the missing pqueries)
 * `chunks/`		based on the chunk placement files from 04_epa, but with the broken pqueries removed (which then become the missing ones)
 * `missing/`	list of missing sequences per sample. should be (almost) empty.
 * `samples/`	the actual sample jplace files, created using a backmapping from the map files, using the chunks jplace files as source

 * `check.uniq`	condesed output of the check/*.diff files, which lists all sequence names that are in the maps, but not in the samples.
 * `chunknames.csv`	lists all chunks and their original mapping name. those were lost during papara and pepa steps.
 * `missing.uniq`	list of all sequence names that were missing in the chunks while creating the samples. those should be the ones removed in the cleaning step
 * `overview.csv`	list of samples and their number of pqueries

The result are placement files for all needed combinations of trees and datasets:

 * `bv_Bact_Constr_backbone`
 * `bv_Bact_Unconstr_backbone`
 * `bv_Gen_Constr_backbone`
 * `bv_Gen_Unconstr_backbone`
 * `hmp_Bact_Constr_backbone`
 * `hmp_Bact_Unconstr_backbone`
 * `hmp_Gen_Constr_backbone`
 * `hmp_Gen_Unconstr_backbone`
 * `tara_Euks_Constr_backbone`
 * `tara_Euks_Unconstr_backbone`
 * `tara_Gen_Constr_backbone`
 * `tara_Gen_Unconstr_backbone`
