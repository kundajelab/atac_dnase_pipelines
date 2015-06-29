import argparse
from multiprocessing import Pool

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam
from pybedtools import BedTool


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create TSS enrichment plots for ATAC-seq data.\n'
                     'Uses cuts i.e. Tn5 adjusted read ends'))
    parser.add_argument('read_bam', help='Sorted read BAM')
    parser.add_argument('tss_file', help='BED/GFF file containing TSS regions')
    parser.add_argument('window_halfwidth', help='Bases to extend around TSS',
                        type=int)
    parser.add_argument('output_prefix', help='Output prefix')
    parser.add_argument('-p', '--processes',
                        help=('Number of processes to run in parallel '
                              '(default: 1)'),
                        type=int, default=1)
    parser.add_argument('-q', '--mapq_threshold',
                        help='Minimum MAPQ threshold',
                        type=int,
                        default=20)
    return parser.parse_args()

TN5_OFFSET_PLUS = 4
TN5_OFFSET_MINUS = -5


if __name__ == '__main__':
    args = parse_args()

    aggregate = np.zeros(2*args.window_halfwidth + 1)

    tss_list = BedTool(args.tss_file)
    bam = pysam.Samfile(args.read_bam, 'rb')

    for tss in tss_list:
        window_start = max(tss.start - args.window_halfwidth - TN5_OFFSET_PLUS,
                           0)
        window_end = tss.end + args.window_halfwidth - TN5_OFFSET_MINUS
        for read in bam.fetch(tss.chrom, window_start, window_end):
            if read.mapq < args.mapq_threshold:
                continue

            # Get Tn5 adjusted 5' ends
            if read.is_reverse:
                cut_pos = read.reference_end - 1 + TN5_OFFSET_MINUS
            else:
                cut_pos = read.reference_start + TN5_OFFSET_PLUS

            try:
                aggregate[cut_pos - window_start] += 1
            except IndexError:
                pass

    fig = plt.figure(figsize=(8.0, 5.0))
    plot_x = range(-args.window_halfwidth, args.window_halfwidth + 1)

    background = np.mean(aggregate[1:200])
    enrichments = aggregate / background
    # use a rectangular 50bp kernel for smoothing
    smoothed_enrichments = (np.convolve(aggregate, np.ones(50), 'same') / 50
                            / background)

    plt.plot(plot_x, enrichments, 'k.')
    plt.plot(plot_x, smoothed_enrichments, 'r')
    plt.xlabel('Position relative to center')
    plt.ylabel('Cut enrichment')
    fig.savefig(args.output_prefix + '.png')
    np.savetxt(args.output_prefix + '.values.txt', aggregate, fmt='%d')
    np.savetxt(args.output_prefix + '.enrichments.txt', enrichments, fmt='%.3f')
    plt.close(fig)
