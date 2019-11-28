#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import click

from conf.plot_settings import parse_settings
from utils.reading_input import index_gtf
from utils.reading_input import read_reads_depth_from_count_table
from utils.reading_input import read_transcripts
from utils.sashimi_plot_utils import draw_sashimi_plot
from utils.utils import *


__dir__ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@click.command()
@click.option(
    "-e",
    "--event",
    type=click.STRING,
    required=True,
    help="Event range eg: chr1:100-200:+"
)
@click.option(
    "-i",
    "--input",
    type=click.Path(exists=True),
    required=True,
    help="Path to junctions count table"
)
@click.option(
    "--input-list",
    type=click.Path(exists=True),
    help="""
    Path to tab separated list file\b
    1. the column use to plot sashimi, identical with count table column names
    2. optional, the alias of 1st column
    3. additional columns
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
    help="Path to output graph file",
    show_default=True
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
    "-d",
    "--dpi",
    default=300,
    type=click.IntRange(min=1, clamp=True),
    help="The resolution of output file",
    show_default=True
)
@click.option(
    "--indicator-lines",
    default=None,
    type=click.STRING,
    help="Where to plot additional indicator lines, comma separated int"
)
@click.option(
    "--share-y",
    default=False,
    is_flag=True,
    type=click.BOOL,
    help="Whether different sashimi plots shared same y axis"
)
@click.option(
    "--no-gene",
    is_flag=True,
    type=click.BOOL,
    help="Do not show gene id next to transcript id"
)
@click.option(
    "--color-factor",
    default=1,
    type=click.IntRange(min=1),
    help="Index of column with color levels (1-based)",
    show_default=True
)
def no_bam(
        event,
        input,
        input_list,
        gtf,
        output,
        config,
        threshold,
        indicator_lines,
        share_y,
        no_gene,
        color_factor,
        dpi
):
    u"""
    This function is used to plot sashimi without BAM file
    \f
    :param event:
    :param input:
    :param input_list:
    :param gtf: path to gtf file
    :param output: path to output file
    :param event: event id, chr:100-200-100-200:+ etc
    :param config: path to config file, default using settings.ini file under this suite of scripts
    :param threshold: filter out low abundance junctions
    :param indicator_lines: draw vertical lines in sashimi to show the spliced sites
    :param share_y: make different plots use same y axis boundary
    :param no_gene: do not show gene id
    :param color_factor: 1-based index, only work with bam list
    :param dpi: output file resolution
    :return:
    """
    sashimi_plot_settings = parse_settings(config)

    colors = sashimi_plot_settings["colors"]
    color_index = 0
    tmp_color = {}
    required_cols = {}
    if input_list:
        with open(input_list) as r:
            for line in r:
                lines = re.split(r"\t| {2,}", line.strip())

                if not lines:
                    continue

                if len(lines) > 1:
                    required_cols[lines[0]] = lines[1]
                else:
                    required_cols[lines[0]] = clean_star_filename(lines[0])

                if lines[color_factor - 1] not in tmp_color.keys():
                    tmp_color[required_cols[lines[0]]] = colors[color_index % len(colors)]
                    color_index += 1

    splice_region = get_sites_from_splice_id(event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        region=splice_region
    )

    reads_depth = read_reads_depth_from_count_table(
        count_table=input,
        splice_region=splice_region,
        required=required_cols,
        threshold=threshold,
        colors=tmp_color
    )

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region,
        no_bam=True,
        show_gene=not no_gene,
        dpi=dpi
    )
