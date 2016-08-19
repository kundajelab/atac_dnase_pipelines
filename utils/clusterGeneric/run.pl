#!/usr/bin/perl

use POSIX;

die "Error: Missing arguments.\nUsage: run.pl timeout cpus mem queue saveStdout saveStderr cmd arg1 ... argN\n" if $#ARGV < 6 ;

$timeout = shift @ARGV;
$cpus = shift @ARGV;
$mem = shift @ARGV;
$queue = shift @ARGV;
$saveStdout = shift @ARGV;
$saveStderr = shift @ARGV;
$cmd = join(' ', @ARGV);

$qsub = "sbatch --export=ALL ";
$qsub .= "-n 1 --ntasks-per-node=1 --cpus-per-task=$cpus " if( $cpus > 0 );
if( $mem > 0 ) {
	$mem = ceil($mem/1000000); # MB
	$qsub .= "--mem-per-cpu $mem ";
}
if( $timeout > 0 ) {
	$timeout = ceil($timeout/60); # minute
	$qsub .= "-t $timeout ";
}

$pid = open QSUB, " | $qsub";
die "Cannot run command '$qsub'\n" if ! kill(0, $pid); # Check that process exists
print QSUB "#!/bin/sh \n";	# SLURM sbatch needs this shebang...
print QSUB "$cmd\n";		# Send cluster's task via qsub's STDIN
close QSUB;

exit(0);

