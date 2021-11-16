#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
from multiprocessing import cpu_count

import click
from conf.plot_settings import parse_settings
from ioutils.reading_input import (index_gtf, read_reads_depth_from_bam,
                                   read_reads_depth_from_count_table,
                                   read_transcripts)
from ioutils.utils import *

from plot.sashimi_plot_utils import draw_sashimi_plot

__dir__ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def darw_without_bam(
    bam,
    count_table,
    splice_region,
    colors,
    color_factor,
    threshold,
    output,
    sashimi_plot_settings,
    no_gene,
    dpi,
    distance_ratio,
    title
):
    color_index = 0
    tmp_color = {}
    required_cols = {}
    if bam:
        with open(bam) as r:
            for line in r:
                lines = re.split(r"\t| {2,}", line.strip())

                if not lines:
                    continue

                if len(lines) > 1:
                    required_cols[lines[0]] = lines[1]
                else:
                    required_cols[lines[0]] = clean_star_filename(lines[0])

                if lines[color_factor - 1] not in tmp_color.keys():
                    tmp_color[required_cols[lines[0]]
                              ] = colors[color_index % len(colors)]
                    color_index += 1

    reads_depth = read_reads_depth_from_count_table(
        count_table=count_table,
        splice_region=splice_region,
        required=colors,
        threshold=threshold,
        colors=colors
    )

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region,
        no_bam=True,
        show_gene=not no_gene,
        dpi=dpi,
        distance_ratio=distance_ratio,
        title=title,
    )


def draw_default(
    bam, splice_region,
    color_factor, colors,
    share_y, share_y_by,
    barcode,
    sort_by_color,
    threshold, threshold_of_reads,
    log, process, reads,
    barcode_tag, save_depth,
    output,
    customized_junction,
    sashimi_plot_settings,
    no_gene, dpi,
    distance_ratio,
    title, show_side,
    stack, side_strand,
    strandness,
    is_atac
):
    bam_list, shared_y = prepare_bam_list(
        bam, color_factor, colors, 
        share_y_by, barcodes=barcode,
        is_atac = is_atac
    )

    if sort_by_color:
        bam_list = sorted(bam_list, key=lambda x: x.color)

    reads_depth = read_reads_depth_from_bam(
        bam_list=bam_list,
        splice_region=splice_region.copy(),
        threshold=threshold,
        threshold_of_reads=threshold_of_reads,
        log=log,
        n_jobs=process,
        reads=reads,
        barcode_tag=barcode_tag,
        strandness=strandness,
        is_atac = is_atac,
        stack = stack
    )

    if save_depth:
        with open(output, "w+") as w:
            for key, value in reads_depth.items():
                for val in value:
                    w.write("{},{}\n".format(key.to_csv(), val))
        exit(0)

    # set shared y
    if share_y:
        assign_max_y(shared_y.values(), reads_depth)

    # read customized junctions
    if customized_junction and os.path.exists(customized_junction):
        customized_junction = read_reads_depth_from_count_table(
            customized_junction,
            splice_region=splice_region,
            required=None,
            colors=colors
        )

        temp_customized_junction = {
            k.alias: v for k, v in customized_junction.items()}

        for key, value in reads_depth.items():
            temp = temp_customized_junction.get(
                key.alias,
                customized_junction.get(os.path.basename(key.path), None)
            )

            if temp:
                value.add_customized_junctions(temp)

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region.copy(),
        no_bam=False,
        show_gene=not no_gene,
        dpi=dpi,
        logtrans=log,
        distance_ratio=distance_ratio,
        title=title,
        stack=stack,
        show_side=show_side,
        side_strand_choice={"+": "plus", "-": "minus"}.get(side_strand)
    )


@click.command()
@click.option(
    "-e",
    "--event",
    type=click.STRING,
    required=True,
    help="Event range eg: chr1:100-200:+"
)
@click.option(
    "-b",
    "--bam",
    type=click.Path(exists=True),
    required=True,
    help="""
    Path to input BAM file. \b

    Or a tab separated text file, \b 
    - first column is path to BAM file, \b
    - second column is BAM file alias(optional) \b
    \b
    Path to tab separated list file\b
    1. the column use to plot sashimi, identical with count table column names \b
    2. optional, the alias of 1st column \b
    3. additional columns \b
    """
)
@click.option(
    "-c",
    "--count-table",
    type=click.Path(),
    help="""
    Path to input count table file. \b
    To make sashimi plot without bam
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
    "-T",
    "--threshold-of-reads",
    default=0,
    type=click.IntRange(min=0, clamp=True),
    help="Threshold to filter low abundance reads for stacked plot",
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
    help="""
    Where to plot additional indicator lines, comma separated int, sites occured multiple times will highlight in red \b
    Or \b
    Path to file contains indicator lines, \b
    1st column is the line site
    2nd column is transcript id
    3rd column is the weights
    """
)
@click.option(
    "--share-y",
    default=False,
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether different sashimi plots shared same y axis
    """
)
@click.option(
    "--no-gene",
    is_flag=True,
    type=click.BOOL,
    help="Do not show gene id next to transcript id"
)
@click.option(
    "--color-factor",
    type=str,
    help="""
    The index of specific column in --bam or path to color settings, 2 columns are required, first if key of bam or cell group, second column is color
    """,
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
    type=click.IntRange(min=1, max=cpu_count()),
    default=1,
    help="""
    How many cpu to use \b
    """
)
@click.option(
    "-f",
    "--genome",
    type=click.Path(),
    default=None,
    help="""
    Path to genome fasta \b
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
@click.option(
    "--stack",
    default=False,
    is_flag=True,
    type=click.BOOL,
    help="Whether to draw stacked reads"
)
@click.option(
    "--share-y-by",
    type=click.INT,
    default=-1,
    help="""
    Index of column with share y axis (1-based), Need --share-y. 
    For example, first 3 bam files use same y axis, and the rest use another
    """,
    show_default=True
)
@click.option(
    "--remove-empty-gene",
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether to plot empty transcript \b
    """
)
@click.option(
    "--distance-ratio",
    type=click.FLOAT,
    default=0.3,
    help="distance between transcript label and transcript line",
    show_default=True
)
@click.option(
    "--title",
    type=click.STRING,
    default=None,
    help="Title",
    show_default=True
)
@click.option(
    "--save-depth",
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether to save reads depth to file, 
    the last 3 columns are chrom, position and depth,
    The same pos will repeated multiple times for joyplot in R
     \b
    """
)
@click.option(
    "--barcode",
    type=click.Path(),
    help="""
    Path to barcode list file, 
    At list  three columns were required,
    1st The name of bam file; 2nd the barcode;
    3rd The group label
     \b
    """
)
@click.option(
    "--barcode-tag",
    type=click.STRING,
    default="CB",
    help="""
    The default cell barcode tag label
     \b
    """
)
@click.option(
    "--reads",
    type=click.Choice(["All", "R1", "R2"]),
    default="All",
    help="""
    Whether filter R1 or R2
     \b
    """
)
@click.option(
    "--show-side",
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether to draw additional side plot, 
     \b
    """
)
@click.option(
    "--side-strand",
    type=click.Choice(["All", "+", "-"]),
    default="All",
    help="""
    which strand kept for side plot, default use all
     \b
    """
)
@click.option(
    "--show-id",
    is_flag=True,
    help="""
    which show gene id or gene name
    \b
    """
)
@click.option(
    "-S", "--strand-specific",
    is_flag=True,
    help="""
    only show transcripts and reads of input region
    \b
    """
)
@click.option(
    "--sc-atac",
    is_flag=True,
    type=click.BOOL,
    help="""
    Whether input file is fragments of scATAC-seq
     \b
    """
)
def plot(
        bam: str, event: str, gtf: str, output: str,
        config: str, threshold: int, indicator_lines: str,
        share_y: bool, no_gene: bool, color_factor: str,
        dpi: int, log: str, customized_junction: str,
        process: int, sort_by_color: bool, share_y_by: int,
        remove_empty_gene: bool, distance_ratio: float,
        title: str, genome: str, save_depth: str, stack: bool,
        threshold_of_reads: int, barcode: str, barcode_tag: str,
        reads: str, show_side: bool, side_strand: str, count_table: str, strand_specific: bool,
        show_id: bool = False, sc_atac: bool = False
):
    u"""
    This function is used to plot single sashimi plotting
    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :param ctx: passed parameters from main
    :param bam: list of input BAM files
    :param event: event id, chr:100-200-100-200:+ etc
    :param bam: list of input BAM files
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
    :param log: whether to perform y axis log transform
    :param customized_junction: add customized junction to plot
    :param process:
    :param sort_by_color:
    :param share_y_by:
    :param remove_empty_gene:
    :param distance_ratio:
    :param title
    :return:
    """
    try:
        log = int(log)
    except ValueError:
        pass

    out_dir = os.path.dirname(os.path.abspath(output))

    # check reads to kept
    reads = None if reads == "All" else reads == "R1"

    try:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    except IOError as err:
        print(err)
        print("Create output directory failed, please check %s" % out_dir)
        exit(err)

    sashimi_plot_settings = parse_settings(config)

    colors = sashimi_plot_settings["colors"]

    splice_region = get_sites_from_splice_id(
        event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        genome=genome,
        region=splice_region,
        show_id=show_id,
        strandness=not strand_specific,
    )

    if remove_empty_gene:
        splice_region.remove_empty_transcripts()

    if count_table:
        darw_without_bam(
            bam,
            count_table,
            splice_region,
            colors,
            color_factor,
            threshold,
            output,
            sashimi_plot_settings,
            no_gene,
            dpi,
            distance_ratio,
            title
        )
    else:
        draw_default(
            bam, splice_region,
            color_factor, colors,
            share_y, share_y_by,
            barcode,
            sort_by_color,
            threshold, threshold_of_reads,
            log, process, reads,
            barcode_tag, save_depth,
            output,
            customized_junction,
            sashimi_plot_settings,
            no_gene, dpi,
            distance_ratio,
            title, show_side,
            stack, side_strand,
            not strand_specific,
            is_atac = sc_atac
        )


if __name__ == '__main__':
    pass
