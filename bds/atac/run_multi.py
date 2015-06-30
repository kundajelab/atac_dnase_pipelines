#!/usr/bin/env python
import argparse
import errno
import os
import sys

from sh import cd
from sh import bds


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            sys.stderr.write('Warning: directory %s exists.\n' % path)
            pass
        else:
            raise


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Run ATAC pipeline on multiple samples.'))
    parser.add_argument('sample_file', help='TSV file containing sample info')
    parser.add_argument('output_directory', help='Output directory')
    parser.add_argument('-s', '--system',
                        help='BDS system parameter (e.g. local/cluster)')
    parser.add_argument('-t', '--threads',
                        help='Number of threads (default: 1)',
                        type=int, default=1)
    return parser.parse_args()


ATAC_BDS = '/users/leepc12/code/pipelines/bds/atac/atac.bds'
MOD_DEF = ('"bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; '
           'ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; '
           'texlive/2013"')
GENOME2GENOME_SIZE = {
    'hg19': 'hs',
    'mm9': 'mm',
}
GENOME2CHROM_SIZE = {
    'hg19': '/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes',
    'mm9': '/srv/scratch/leepc12/mm9/mm9.chrom.sizes',
}


if __name__ == "__main__":
    args = parse_args()

    with open(args.sample_file) as fp:
        for line in fp:
            line = line.strip().split('\t')
            sample_id, description, read1, read2, genome, index, tssfile = line

            sample_output_dir = '%s-%s' % (sample_id, description)

            cd(args.output_directory)
            mkdir_p(sample_output_dir)
            genome_size = GENOME2GENOME_SIZE[genome]
            chrom_size = GENOME2CHROM_SIZE[genome]
            bds('-s', args.system, ATAC_BDS, index, read1, read2, args.threads,
                genome_size, chrom_size, tssfile, sample_output_dir,
                '-mod', MOD_DEF)
