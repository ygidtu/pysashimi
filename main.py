#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import os
import re

import click

from collections import namedtuple

from src.reading_input import SpliceRegion
from src.plot_settings import parse_settings
from src.reading_input import read_reads_depth, read_transcripts, index_gtf, is_bam
from src.sashimi_plot_utils import draw_sashimi_plot

__dir__ = os.path.dirname(__file__)
VERSION = "1.0.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


bam_info = namedtuple("bam_info", ["alias", "title", "label", "path"])


def get_sites_from_splice_id(string, sep="@", span=0, indicator_lines=None):
    u"""
    get splice range from splice id
    :param string: splice id
    :param sep: the separator between different junctions
    :param span: the range to span the splice region, int -> bp; float -> percentage
    :param indicator_lines: bool
    :return: chromosome, start, end, strand
    """
    u"""
     get splice range from splice id
     :param string: splice id
     :return: chromosome, start, end, strand
     """
    if re.search(r"[\w\.]:(\d+-?){2,}:[+-]", string):
        chromosome, sites, strand = string.split(":")
    elif re.search(r"[\w\.]:(\d+-?){2,}[+-]", string):
        chromosome, sites = string.split(":")
        sites, strand = sites[:-1], sites[-1]
    else:
        chromosome, sites = string.split(":")
        strand = "."


    sites = sites.split("-")

    return SpliceRegion(
        chromosome=chromosome,
        start=sites[0],
        end=sites[-1],
        strand=strand,
        events=string,
        sites=[int(x) for x in indicator_lines.split(",")] if indicator_lines else None
    )


@click.command(
    "Plot sashimi plot",
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "-e"
    "--event",
    type=click.STRING,
    required=True,
    help="Event range eg: chr1:100-200:+"
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
@click.option(
    "-t",
    "--threshold",
    default=0,
    type=click.IntRange(min=0, clamp=True),
    help="Threshold to filter low abundance junctions",
    show_default=True
)
@click.option(
    "--indicator-lines",
    default=None,
    type=click.STRING,
    help="Where to plot additional indicator lines, comma separated int"
)
@click.version_option(VERSION, message="Current version %(version)s")
def main(
        bam,
        gtf,
        output,
        event,
        config,
        threshold,
        indicator_lines
):
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
    :param threshold:
    :param indicator_lines:
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

    clean_bam_filename = lambda x: re.sub("[_.]Aligned.sortedByCoord.out.bam", "", os.path.basename(x))

    if is_bam(bam):
        bam_list = [
            bam_info(
                path=bam,
                alias=clean_bam_filename(bam),
                title=None,
                label=None
            )
        ]

    else:
        bam_list = []
        with open(bam) as r:
            for line in r:
                lines = line.split()

                if len(lines) > 1:
                    tmp = bam_info(
                        path=lines[0],
                        alias=lines[1],
                        title=None,
                        label=None
                    )
                else:
                    tmp = bam_info(
                        path=bam,
                        alias=clean_bam_filename(bam),
                        title=None,
                        label=None
                    )
                bam_list.append(tmp)

    sashimi_plot_settings = parse_settings(config)

    splice_region = get_sites_from_splice_id(event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        region=splice_region
    )

    reads_depth = read_reads_depth(
        bam_list=bam_list,
        splice_region=splice_region,
        threshold=threshold
    )

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region
    )


if __name__ == '__main__':
    main()
    pass
