#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=24:00:00                  # Request the time for the job
#
#@ job_type = parallel
#@ class = micro
#@ node = 1
#@ total_tasks = 28
#@ island_count = 1
#@ energy_policy_tag = raxml_nt_gamma
#@ minimize_time_to_solution = yes
## other example
#@ wall_clock_limit = 48:00:00
#@ job_name = 01_trees_Gen_Unconstr_backbone
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = 01_trees/Gen_Unconstr_backbone/sREPLACE_SEED
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=lucas.czech@h-its.org
#@ queue

SEED=REPLACE_SEED

BASEDIR=/path/to/here
NUM_TASKS=28

ALI=${BASEDIR}/11_single_seqs/00_reference/Archaea_sequences.cleaned.fasta.reduced

RAXML=${BASEDIR}/software/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3

${RAXML} -f o -p ${SEED} -m GTRGAMMAX -s ${ALI} -n s${SEED} -T ${NUM_TASKS}
