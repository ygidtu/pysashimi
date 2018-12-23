#!/usr/bin/env python3
#-*-coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Tests
"""
from tqdm import tqdm
import multiprocessing as mp
import os
import sys
import argparse as ap
import pickle
# import matplotlib.pyplot as plt
# from matplotlib.path import Path
# from matplotlib.patches import PathPatch
# from pylab import linspace

from src.sashimi_plot_utils import draw_sashimi_plot
from src.reading_input import read_reads_depth, read_transcripts, index_gtf, logger
from src.plot_settings import parse_settings


__dir__ = os.path.dirname(__file__)


'''
def test_path(show_axis, chrom, strand, font_size, nxticks=4):
    pts = [
        (0, 0),
        (0, 1),
        (5, 1),
        (5, 0)
    ]

    a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    p = PathPatch(a, ec="red", lw=1, fc='none')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(p)
    ax.set_xlim(-2, 7)
    ax.set_ylim(-1, 5)
    ax.text(2.5, 1.1, "10")

    # 去除一半，右侧和上侧的坐标轴
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    # if show_axis:
    #     ax.xaxis.set_ticks_position('bottom')
    #     ax.xlabel('Genomic coordinate (%s), "%s" strand'%(chrom,
    #                                                    strand),
    #            fontsize=font_size)
    #     max_graphcoords = max(graphcoords) - 1
    #     ax.xticks(
    #         linspace(0, max_graphcoords, nxticks),
    #            [graphToGene[int(x)] for x in \
    #             linspace(0, max_graphcoords, nxticks)],
    #            fontsize=font_size
    #     )
    # else:
    #     ax.spines['bottom'].set_color('none')
    #     ax.xticks([])
    #
    # ax.xlim(0, max(graphcoords))

    plt.show()
'''


def get_sites_from_splice_id(string):
    u"""
    get splice range from splice id
    :param string: splice id
    :return: chromosome, start, end, strand
    """
    chromosome, sites, strand = string.split(":")
    sites = sites.split("-")
    return chromosome, int(sites[0]), int(sites[-1]), strand


def command_line_args():
    u"""
    generate command line parameters
    :return:
    """
    parser = ap.ArgumentParser(description="Sashimi plot")

    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        # required=True,
        help="Path to bam files"
    )
    parser.add_argument(
        "-g",
        "--gtf",
        type=str,
        help="Path to gtf file"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to output file"
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default=os.path.join(__dir__, "new_src/settings.ini"),
        # required=True,
        help="Path to setting file"
    )
    parser.add_argument(
        "-e",
        "--event",
        type=str,
        required=True,
        help="splice event id"
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    else:
        try:
            args = parser.parse_args(sys.argv[1:])
            print(vars(args))
            return args
        except ap.ArgumentError as err:
            print(err)
            parser.print_usage()
            exit(err)


def test_sashimi(args):
    u"""
    Created by Zhang yiming at 2018.12.19
    This function is used to test the function of sashimi plotting
    :param bam: list of input BAM files
    :param gtf: path to gtf file
    :param output: path to output file
    :param event: event id, chr:100-200-100-200:+ etc
    :param config: path to config file, default using settings.ini file under this suite of scripts
    :return:
    """
    # bam, gtf, output, event, config = args
    bam, gtf, output, event, config = [
            "/Volumes/WD/tests/splice_plot/L02_102_chr1.bam",
            "/Volumes/WD/tests/splice_plot/Homo_sapiens.GRCh38.93.sorted.gtf",
            "/Volumes/WD/tests/splice_plot/test.pdf",
            "1:8062-29320:-",
            os.path.join(__dir__, "src/settings.ini")
    ]

    output = os.path.abspath(output)
    out_dir = os.path.dirname(output)

    try:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    except IOError as err:
        print("Create output directory failed, please check %s" % out_dir)

    bam_list = []
    for i in bam.split(","):
        if not os.path.exists(i):
            raise ValueError("%s not found" % i)

        bam_list.append(os.path.abspath(i))

    logger.info("check")
    chromosome, start, end, strand = get_sites_from_splice_id(event)

    logger.info("splie")
    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand
    )

    reads_depth = read_reads_depth(
        bam_list=bam_list,
        splice_region=splice_region
        # strand=strand
    )

    """
    output_file_path,
    settings,
    var_pos,
    average_depths_dict,
    mRNAs_object,
    """
    sashimi_plot_settings = parse_settings(config)

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region
    )


def test_in_batch(infile, output, gtf, bam, config, n_job):
    u"""
    test plot in batch
    :return:
    """
    args = []
    with open(infile) as r:
        for line in r:
            if line.startswith("#"):
                continue
            lines = line.split()

            tmp = [
                bam,
                gtf,
                os.path.join(output, lines[3] + ".pdf"),
                lines[3],
                config,
            ]

            # if os.path.exists(tmp[2]):
            #     continue

            args.append(tmp)

    with mp.Pool(processes=n_job) as pool:
        list(tqdm(pool.imap(test_sashimi, args), total=len(args)))


if __name__ == '__main__':
    # main()

    # with open(os.path.join(__dir__, "test_files/test.p"), "rb") as r:
    #     data = pickle.load(r)
    #
    # print(data)

    # args = command_line_args()

    # ,/Volumes/WD/tests/splice_plot/L02_99_chr1.bam
    logger.info("start")
    test_sashimi([
            "/Volumes/WD/tests/splice_plot/L02_102_chr1.bam",
            "/Volumes/WD/tests/splice_plot/Homo_sapiens.GRCh38.93.sorted.gtf",
            "/Volumes/WD/tests/splice_plot/test.pdf",
            "1:8062-29320:-",
            os.path.join(__dir__, "src/settings.ini")
    ])

    # test_sashimi(
    #     [
    #         args.bam,
    #         args.gtf,
    #         args.output,
    #         args.event,
    #         args.config,
    #     ]
    # )

    # test_in_batch(
    #     infile=args.event,
    #     config=args.config,
    #     gtf=args.gtf,
    #     bam=args.bam,
    #     output=args.output,
    #     n_job=12
    # )


