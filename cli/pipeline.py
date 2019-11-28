#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import click

from multiprocessing import cpu_count

from conf.plot_settings import parse_settings
from utils.reading_input import index_gtf
from utils.reading_input import read_reads_depth_from_bam
from utils.reading_input import read_reads_depth_from_count_table
from utils.reading_input import read_transcripts
from utils.sashimi_plot_utils import draw_sashimi_plot
from utils.utils import *


__dir__ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@click.command()
@click.option(
    "-i",
    "--input",
    type=click.Path(exists=True),
    required=True,
    help="Path to the meta info [xlsx]"
)
@click.option(
    "-s",
    "--span",
    default="100",
    type=click.STRING,
    help="""
        To span the input region, \b
        int -> span corresponding bp \b
        float -> span by percentage of input region
    """,
    show_default=True
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
    help="Path to output directory",
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
    is_flag=True,
    default=False,
    type=click.BOOL,
    help="Where to plot additional indicator lines"
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
@click.option(
    '--log',
    type=click.Choice(["0", "2", "10", "zscore"]),
    default="0",
    help="y axis log transformed, 0 -> not log transform; 2 -> log2; 10 -> log10"
)
@click.option(
    "--customized-junction",
    type=click.STRING,
    default=None,
    help="""
    Path to junction table column name needs to be bam name or bam alias. \b
    """
)
@click.option(
    "-p",
    "--process",
    type=click.IntRange(min=1, max = cpu_count()),
    default=1,
    help="""
    How many cpu to use \b
    """
)
@click.option(
    "--sort-by-color",
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether sort input bam order, for better looking \b
    """
)
def pipeline(
        input,
        span,
        gtf,
        output,
        config,
        threshold,
        indicator_lines,
        share_y,
        no_gene,
        color_factor,
        dpi,
        log,
        customized_junction,
        process,
        sort_by_color
):
    u"""

    This function is used to plot sashimi based on specific meta info

    required a specific format of input file

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :param input: input file in specific format
    :param span: str, but must be int or float
    :param gtf: path to gtf file
    :param output: path to output file
    :param config: path to config file, default using settings.ini file under this suite of scripts
    :param threshold: filter out low abundance junctions
    :param indicator_lines: draw vertical lines in sashimi to show the spliced sites
    :param share_y: make different plots use same y axis boundary
    :param no_gene: do not show gene id
    :param color_factor: 1-based index, only work with bam list
    :param dpi: output file resolution
    :param log: whether to perform y axis log transform
    :param customized_junction: add customized junction to plot
    :param process:
    :param sort_by_color
    :return:
    """

    try:
        log = int(log)
    except ValueError:
        pass

    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except IOError as err:
        print(err)
        print("Create output directory failed, please check %s" % output)
        exit(err)

    sashimi_plot_settings = parse_settings(config)

    data, bam_list = read_info_from_xlsx(
        xlsx=input,
        color_factor=color_factor,
        colors=sashimi_plot_settings["colors"]
    )

    if sort_by_color:
        bam_list = sorted(bam_list, key=lambda x: x.color)

    coords = get_merged_event(data.keys(), span=span, indicator_lines=indicator_lines)

    for k, v in coords.items():

        # v = assign_events(v)

        # for merged, separate in v.items():
        for region in v:

            splice_region = read_transcripts(
                gtf_file=index_gtf(input_gtf=gtf),
                region=region.copy()
            )

            reads_depth = read_reads_depth_from_bam(
                bam_list=bam_list,
                splice_region=splice_region.copy(),
                threshold=threshold,
                log=log,
                n_jobs=process
            )

            # read customized junctions
            if customized_junction and os.path.exists(customized_junction):
                customized_junction = read_reads_depth_from_count_table(
                    customized_junction,
                    splice_region=splice_region,
                    required=None,
                    colors=sashimi_plot_settings["colors"]
                )

                for key, value in reads_depth.items():
                    temp = customized_junction.get(
                        key.alias,
                        customized_junction.get(os.path.basename(key.path), None)
                    )

                    if temp:
                        value.add_customized_junctions(temp)

            # for sep in separate:

            tmp_reads_depth_dict = {}

            # add label to read_depth
            for i, j in reads_depth.items():
                tmp_reads_depth_of_bam = j.get_read_depth(region)

                i = i._replace(label=data[region.ori][i.title])
                tmp_reads_depth_dict[i] = tmp_reads_depth_of_bam

            draw_sashimi_plot(
                output_file_path=os.path.join(output, region.events + ".pdf"),
                settings=sashimi_plot_settings,
                average_depths_dict=tmp_reads_depth_dict,
                splice_region=splice_region.get_region(region),
                no_bam=False,
                show_gene=not no_gene,
                dpi=dpi,
                log=log
            )

