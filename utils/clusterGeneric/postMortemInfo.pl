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
# The following command is executed in order to get information of a recently 
# finished jobId. This information is typically used for debuging and it added
# to bds's output.
#
# Script's output: 
#     The output is not parsed, it is stored and later shown 
#     in bds's report. Is should contain information relevant 
#     to the job's execution (e.g. "qstat -f $jobId" or 
#     "checkjob -v $jobId")
#
# Command line arguments: 
#     jobId: This is the jobId returned as the first line in 'clusterGenericRun' 
#           script (i.e. the jobID provided by the cluster management system)
#
#                                                                Pablo Cingolani
#-------------------------------------------------------------------------------

#---
# Parse command line arguments
#---
die "Error: Missing arguments.\nUsage: postMortemInfo.pl jobId\n" if $#ARGV < 0 ;
$jobId = shift @ARGV;

#---
# Execute cluster command to show task details
#---
$exitCode = system "squeue -j $jobId";

# OK
exit($exitCode);

