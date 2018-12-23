#!/usr/bin/env python3
#-*-coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import re
import click
import multiprocessing as mp
import os
import sys
import argparse as ap

from src.sashimi_plot_utils import draw_sashimi_plot
from src.reading_input import read_reads_depth, read_transcripts, index_gtf, is_bam
from src.plot_settings import parse_settings


__dir__ = os.path.dirname(__file__)
VERSION = "1.0.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def get_sites_from_splice_id(string):
    u"""
    get splice range from splice id
    :param string: splice id
    :return: chromosome, start, end, strand
    """
    chromosome, sites, strand = string.split(":")
    sites = sites.split("-")
    return chromosome, int(sites[0]), int(sites[-1]), strand


@click.command(
    "Plot sashimi plot",
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "-b",
    "--bam",
    type=click.Path(exists=True),
    help="""
    Path to input BAM file. \b
    
    Or a tab separated text file, \b 
    - first column is path to BAM file, \b
    - second column is BAM file alias(optional)
    """
)
@click.option(
    "-g",
    "--gtf",
    type=click.Path(exists=True),
    help="Path to gtf file, both transcript and exon tags are necessary"
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    help="Path to output graph file"
)
@click.option(
    "--config",
    default=os.path.join(__dir__, "settings.ini"),
    type=click.Path(),
    help="Path to config file, contains graph settings of sashimi plot",
    show_default=True
)
@click.version_option(VERSION, message="Current version %(version)s")
def main(bam, gtf, output, event, config):
    u"""
    Created by Zhang yiming at 2018.12.19

    \b
    This function is used to test the function of sashimi plotting

    \f
    :param bam: list of input BAM files
    :param gtf: path to gtf file
    :param output: path to output file
    :param event: event id, chr:100-200-100-200:+ etc
    :param config: path to config file, default using settings.ini file under this suite of scripts
    :return:
    """

    output = os.path.abspath(output)
    out_dir = os.path.dirname(output)

    try:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    except IOError as err:
        print(err)
        print("Create output directory failed, please check %s" % out_dir)

    clean_bam_filename = lambda x: re.sub("[_\.]Aligned.sortedByCoord.out.bam", "", os.path.basename(x))

    if is_bam(bam):
        bam_list = {
            clean_bam_filename(bam): os.path.abspath(bam)
        }
    else:
        bam_list = {}
        with open(bam) as r:
            for line in r:
                lines = line.split()

                if len(lines) > 1:
                    bam_list[lines[1]] = lines[0]
                else:
                    bam_list[clean_bam_filename[lines[0]]] = os.path.abspath(lines[0])

    chromosome, start, end, strand = get_sites_from_splice_id(event)

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

    sashimi_plot_settings = parse_settings(config)

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region
    )


if __name__ == '__main__':
    main()
    pass


