#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

cli function to plot sashimi plot
"""
from multiprocessing import cpu_count

import click
import matplotlib as mpl
import matplotlib.font_manager
import numpy as np

from conf.plot_settings import parse_settings
from ioutils.reading_input import (index_gtf, read_transcripts, read_reads_depth_from_bam,
                                   read_reads_depth_from_count_table)
from ioutils.utils import *
from plot.sashimi import save_fig
from src.Bigwig import __CLUSTERING_METHOD__, __DISTANCE_METRIC__
from src.logger import init_logger

np.seterr(all="ignore")

mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42

if any(["Arial" in f.name for f in matplotlib.font_manager.fontManager.ttflist]):
    mpl.rcParams['font.family'] = 'Arial'

__dir__ = os.path.dirname(os.path.abspath(__file__))

VERSION = "1.6.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


init_logger("INFO")


def draw_without_bam(
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
                    tmp_color[required_cols[lines[0]]] = colors[color_index % len(colors)]
                    color_index += 1

    reads_depth = read_reads_depth_from_count_table(
        count_table=count_table,
        splice_region=splice_region,
        required=colors,
        threshold=threshold,
        colors=colors
    )

    save_fig(
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
        strandless,
        sc_atac,
        **kwargs
):
    bam_list, shared_y = prepare_bam_list(
        bam, color_factor, colors, share_y_by,
        barcodes=barcode, is_atac=False,
        show_mean=kwargs.get("show_mean", False)
    )

    # if there is any sc_atac files
    if sc_atac:
        sc_atac_list, shared_y_atac = prepare_bam_list(
            sc_atac, color_factor, colors, share_y_by, barcodes=barcode,
            is_atac=True, show_mean=kwargs.get("show_mean", False)
        )

        bam_list += sc_atac_list
        shared_y.update(shared_y_atac)

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
        strandless=strandless,
        stack=stack,
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

        temp_customized_junction = {k.alias: v for k, v in customized_junction.items()}

        for key, value in reads_depth.items():
            temp = temp_customized_junction.get(
                key.alias,
                customized_junction.get(os.path.basename(key.path), None)
            )

            if temp:
                value.add_customized_junctions(temp)

    save_fig(
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
        side_strand_choice={"+": "plus", "-": "minus"}.get(side_strand),
        bigwigs=kwargs.get("bigwigs")
    )


@click.command(context_settings=CONTEXT_SETTINGS, )
@click.version_option(VERSION, message="Current version %(version)s")
@click.option(
    "-e", "--event", type=click.STRING, required=True, show_default=True,
    help="Event range eg: chr1:100-200:+"
)
@click.option(
    "-b", "--bam", type=click.Path(), required=True, show_default=True, default="",
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
    "--sc-atac", type=click.Path(), show_default=True, default="",
    help="The list of scATAC-seq fragments, same format with --bam"
)
@click.option(
    "--bigwigs", type=click.Path(), show_default=True, default="",
    help="The list of bigWig files for signal display, same format with --bam"
)
@click.option(
    "-c", "--count-table", type=click.Path(), show_default=True, default="",
    help="Path to input count table file. To make sashimi plot without bam"
)
@click.option(
    "-g", "--gtf", type=click.Path(exists=True), show_default=True,
    help="Path to gtf file, both transcript and exon tags are necessary"
)
@click.option("-o", "--output", type=click.Path(), show_default=True, help="Path to output graph file")
@click.option(
    "--config", default=os.path.join(os.path.dirname(__dir__), "settings.ini"),
    type=click.Path(), show_default=True,
    help="Path to config file, contains graph settings of sashimi plot"
)
@click.option(
    "-t", "--threshold", default=0, type=click.IntRange(min=0, clamp=True), show_default=True,
    help="Threshold to filter low abundance junctions"
)
@click.option(
    "-T", "--threshold-of-reads", default=0, type=click.IntRange(min=0, clamp=True), show_default=True,
    help="Threshold to filter low abundance reads for stacked plot"
)
@click.option(
    "-d", "--dpi", default=300, type=click.IntRange(min=1, clamp=True), show_default=True,
    help="The resolution of output file"
)
@click.option(
    "--indicator-lines", default=None, type=click.STRING, show_default=True,
    help="""
    Where to plot additional indicator lines, comma separated int, sites occurred multiple times will highlight in red\b
    Or \b
    Path to file contains indicator lines, \b
    1st column is the line site
    2nd column is transcript id
    3rd column is the weights
    """
)
@click.option(
    "--share-y", default=False, is_flag=True, type=click.BOOL, show_default=True,
    help="Whether different sashimi plots shared same y axis"
)
@click.option(
    "--no-gene", is_flag=True, type=click.BOOL, show_default=True,
    help="Do not show gene id next to transcript id"
)
@click.option(
    "--color-factor", type=str, show_default=True,
    help="""
    The index of specific column in --bam or path to color settings, 
    2 columns are required, first if key of bam or cell group, second column is color
    """
)
@click.option(
    '--log', type=click.Choice(["0", "2", "10", "zscore"]), default="0", show_default=True,
    help="y axis log transformed, 0 -> not log transform; 2 -> log2; 10 -> log10"
)
@click.option(
    "--customized-junction", type=click.STRING, default=None, show_default=True,
    help="Path to junction table column name needs to be bam name or bam alias."
)
@click.option(
    "-p", "--process", type=click.IntRange(min=1, max=cpu_count()), default=1, show_default=True,
    help="How many cpu to use"
)
@click.option("-f", "--genome", type=click.Path(), default=None, show_default=True, help="Path to genome fasta")
@click.option(
    "--sort-by-color", is_flag=True, type=click.BOOL, show_default=True,
    help="Whether sort input bam order, for better looking"
)
@click.option(
    "--stack", default=False, is_flag=True, type=click.BOOL, show_default=True,
    help="Whether to draw stacked reads"
)
@click.option(
    "--share-y-by", type=click.INT, default=-1, show_default=True,
    help="""
    Index of column with share y axis (1-based), Need --share-y. \b
    For example, first 3 bam files use same y axis, and the rest use another
    """
)
@click.option(
    "--remove-empty-gene", is_flag=True, type=click.BOOL, show_default=True,
    help="Whether to plot empty transcript"
)
@click.option(
    "--distance-ratio", type=click.FLOAT, default=0.3, show_default=True,
    help="Distance between transcript label and transcript line"
)
@click.option("--title", type=click.STRING, default=None, help="Title", show_default=True)
@click.option(
    "--save-depth", is_flag=True, type=click.BOOL, show_default=True,
    help="""
    Whether to save reads depth to file, \b
    the last 3 columns are chrom, position and depth,\b
    The same pos will repeated multiple times for joyplot in R
    """
)
@click.option(
    "--barcode", type=click.Path(), show_default=True,
    help="""
    Path to barcode list file, \b
    At list  three columns were required,\b
    1st The name of bam file; 2nd the barcode;\b
    3rd The group label
    """
)
@click.option(
    "--barcode-tag", type=click.STRING, default="CB", show_default=True,
    help="The default cell barcode tag label"
)
@click.option(
    "--reads", type=click.Choice(["All", "R1", "R2"]), default="All", show_default=True,
    help="Whether filter R1 or R2"
)
@click.option(
    "--show-side", is_flag=True, type=click.BOOL, show_default=True,
    help="Whether to draw additional side plot"
)
@click.option(
    "--side-strand", type=click.Choice(["All", "+", "-"]), default="All", show_default=True,
    help="Which strand kept for side plot, default use all"
)
@click.option("--show-id", is_flag=True, show_default=True, help="Which show gene id or gene name")
@click.option(
    "-S", "--strand-specific", is_flag=True, show_default=True,
    help="Only show transcripts and reads of input region"
)
@click.option("--show-mean", is_flag=True, type=click.BOOL, show_default=True, help="Show mean coverage by groups")
@click.option("--focus", type=click.STRING, show_default=True, help="The highlight regions: 100-200:300-400")
@click.option(
    "--stroke", type=click.STRING, show_default=True,
    help="The stroke regions: start1-end1:start2-end2@color-label, draw a stroke line at bottom, default color is red"
)
@click.option("--bigwig-clustering", is_flag=True, show_default=True, help="The clustering the bigwig files")
@click.option(
    "--bigwig-clustering-method", type=click.Choice(__CLUSTERING_METHOD__), default="ward", show_default=True,
    help="The clustering method for the bigwig files"
)
@click.option(
    "--bigwig-distance-metric", type=click.Choice(__DISTANCE_METRIC__), default="euclidean", show_default=True,
    help="The clustering method for the bigwig files"
)
@click.option("--bigwig-scale", is_flag=True, show_default=True, help="The do scale on bigwigs files")
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
        show_id: bool = False, sc_atac: bool = False,
        show_mean: bool = False, focus: str = "",
        bigwigs: str = "",
        bigwig_clustering: bool = False,
        bigwig_clustering_method: str = "ward",
        bigwig_distance_metric: str = "euclidean",
        bigwig_scale: bool = False, stroke: str = ""
):
    u"""
    This function is used to plot single sashimi plotting
    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :return:
    """
    try:
        log = int(log)
    except ValueError:
        pass

    if all([not os.path.exists(x) for x in [bam, sc_atac, bigwigs]]):
        raise ValueError("At least provide --bam, --sc-atac or --bigwigs")

    out_dir = os.path.dirname(os.path.abspath(output))

    # check reads to kept
    reads = None if reads == "All" else reads == "R1"

    try:
        os.makedirs(out_dir, exist_ok=True)
    except IOError as err:
        logger.warning(err)
        logger.error("Create output directory failed, please check %s" % out_dir)
        exit(err)

    sashimi_plot_settings = parse_settings(config)
    colors = sashimi_plot_settings["colors"]

    splice_region = get_sites_from_splice_id(
        event,
        indicator_lines=indicator_lines,
        focus=focus, stroke=stroke
    )

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        genome=genome,
        region=splice_region,
        show_id=show_id,
        strandless=not strand_specific,
    )

    if remove_empty_gene:
        splice_region.remove_empty_transcripts()

    wigs = prepare_bigwig_list(
        bigwigs, splice_region,
        process=process,
        clustering=bigwig_clustering,
        clustering_method=bigwig_clustering_method,
        distance_metric=bigwig_distance_metric,
        do_scale=bigwig_scale
    )

    if count_table:
        draw_without_bam(
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
            sc_atac=sc_atac,
            show_mean=show_mean,
            focus=focus,
            bigwigs=wigs
        )


if __name__ == '__main__':
    pass
