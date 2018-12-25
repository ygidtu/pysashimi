#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import os
import re

import click

from tqdm import tqdm

from src.reading_input import GenomicLoci
from src.plot_settings import parse_settings
from src.reading_input import read_reads_depth, read_transcripts, index_gtf, is_bam
from src.sashimi_plot_utils import draw_sashimi_plot

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
    if re.search(r"[\w\.]:(\d+-?){2,}:[+-]", string):
        chromosome, sites, strand = string.split(":")
    elif re.search(r"[\w\.]:(\d+-?){2,}[+-]", string):
        chromosome, sites = string.split(":")
        sites, strand = sites[:-1], sites[-1]
    else:
        chromosome, sites = string.split(":")
        strand = "."

    sites = sites.split("-")

    return GenomicLoci(
        chromosome=chromosome,
        start=int(sites[0]),
        end=int(sites[-1]),
        strand=strand
    )


def __assign_events__(events):
    u"""
    merge separate events into a huge range
    assign events back to this range
    :param events:
    :return: {merged: [events, events, events]}
    """
    res = {}
    tmp = [events[0]]
    current = events[0]
    for i in events[1:]:
        if current.is_overlap(i):
            current += i
            tmp.append(i)
        else:
            current = i
            res[current] = tmp
            tmp = [i]

    if tmp:
        res[current] = [tmp]

    return res


def plot_batch(event, gtf, bam_list, output, sashimi_plot_settings):
    u"""
    Created by ygidtu at 2018.12.25
    plot multiple events at same time
    :return:
    """

    coords = {}
    for e in event:
        tmp = get_sites_from_splice_id(e)
        tmp_key = "%s#%s" % (tmp.chromosome, tmp.strand)

        tmp_list = coords[tmp_key] if tmp_key in coords.keys() else []
        tmp_list.append(tmp)
        coords[tmp_key] = tmp_list

    for k, v in tqdm(coords.items()):
        tqdm.write(k)

        v = __assign_events__(v)

        for merged, separate in v.items():
            splice_region = read_transcripts(
                gtf_file=index_gtf(input_gtf=gtf),
                chromosome=k.chromosome,
                start=k.start,
                end=k.end,
                strand=k.strand
            )

            reads_depth = read_reads_depth(
                bam_list=bam_list,
                splice_region=splice_region
                # strand=strand
            )

            for sep in separate:
                draw_sashimi_plot(
                    output_file_path=os.path.join(output, str(sep) + ".pdf"),
                    settings=sashimi_plot_settings,
                    average_depths_dict=reads_depth.get_read_depth(sep),
                    splice_region=splice_region.get_region(sep)
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
@click.version_option(VERSION, message="Current version %(version)s")
def main(bam, gtf, output, event, config, threshold):
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
                    bam_list[clean_bam_filename(lines[0])] = os.path.abspath(lines[0])

    genomic = get_sites_from_splice_id(event)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        chromosome=genomic.chromosome,
        start=genomic.start,
        end=genomic.end,
        strand=genomic.strand
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
