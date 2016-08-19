#!/usr/bin/perl

#-------------------------------------------------------------------------------
# BDS generic cluster example
#
# This is a trivial example of the 'cluster generic' interface implementation.
# The commands implemented in this example simply pass the propper arguments 
# to qsub, qdel or qstat commands.
# This is intended as a toy example, since bds can do this directly (but 
# it's a good starting point to extend your own implementation).
#
# This script is executed in order to show the jobID of all jobs currently 
# scheduled in the cluster
#
# Script's output: 
#     This script is expected to print all jobs currently scheduled or 
#     running in the cluster (e.g. qstat). One per line. The FIRST column 
#     should be the jobID (columns are spce or tab separated). Other 
#     columns may exists (but are currently ignored).
#
# Command line arguments: 
#     None
#
#                                                                Pablo Cingolani
#-------------------------------------------------------------------------------

#---
# Execute cluster command to show all tasks
#---
$exitCode = system "squeue";

# OK
exit($exitCode);
