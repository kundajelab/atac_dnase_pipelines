Recursive Tools for BDS reporting/debugging
===================================================================

### General description

The main purpose of all tools here is to copy (or make a symlink of) files on [SRC] to [DEST] recursively. 

The same directory structure as in [SRC] will be created on [DEST]. All dead soft links and empty directory on [DEST] will be automatically removed.

These tools are very useful if you have a web server connected to your computation server and you want to monitor status of your BDS jobs real-time online.


### Examples

```
#Example, copy *.bigwig on [SRC] to [DEST] recursivley.
#The same directory structure as in [SRC] will be created on [DEST].

$recursive_cp.sh [SRC] [DEST] *.bigwig


#Example, make symlinks on [DEST] for *.pdf on [SRC] recursivley.
#The same directory structure as in [SRC] will be created on [DEST].

$recursive_ln.sh [SRC] [DEST] *.pdf


#Example, copy all *.bigwig on [SRC] to [DEST] recursively.
#No directory will be created on [DEST]. 
#Files on [DEST] will have unique names which include directory sturcture as in [SRC].

$recursive_gather_cp.sh [SRC] [DEST] *.bigwig


#Example, make symlinks on [DEST] for *.bam on [SRC]. 
#No directory will be created on [DEST].
#Files on [DEST] will have unique names which include directory sturcture as in [SRC].

$recursive_gather_ln.sh [SRC] [DEST] *.bam
```


### Sync BDS report to web server

It is important to check your BDS job status real-time online. If you have a NFS mounted web directory (eg. /srv/www/kundaje/leepc12/job_status/) connected to your cluster, you can automatically synchronize all files on your working directory to the web directory.
```
./_tools/sync_bds_report_cp.sh [SRC] [DEST]
./_tools/sync_bds_report_ln.sh [SRC] [DEST]
```

You can use a cron job to automate reporting.
```
$crontab -e

# This is an example for nandi cluster. Sync important files every 5 minutes.
# Modify path for sync_bds_report_??.sh, SRC (top of your working directory) and DEST (web directory)

# sync_bds_report_ln.sh generates symlinks on [DEST] for ALL files on [SRC], use this if your web server can see files on [SRC]
*/5 * * * * /users/leepc12/code/pipelines/bds/_tools/sync_bds_report_ln.sh /srv/scratch/leepc12/run /srv/www/kundaje/leepc12/bds_monitor/mitra

# sync_bds_report_cp.sh transfers IMPORTANT files (*.html,*.js,*.log,*.pdf,*.txt) on [SRC] to [DEST], use this if your web server cannot see files on [SRC]
*/5 * * * * /users/leepc12/code/pipelines/bds/_tools/sync_bds_report_cp.sh /srv/scratch/leepc12/run /srv/www/kundaje/leepc12/bds_monitor/nandi
```

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
