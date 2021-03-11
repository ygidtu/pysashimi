
#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import click

from multiprocessing import cpu_count

from conf.plot_settings import parse_settings
from ioutils.reading_input import index_gtf
from ioutils.reading_input import read_reads_depth_from_bam
from ioutils.reading_input import read_transcripts
from ioutils.utils import *
from plot.line_plot_utils import draw_line_plot


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
    "-b",
    "--bam",
    type=click.Path(exists=True),
    help="""
    a tab separated text file, \b 
    - 1st column is path to BAM file, \b
    - 2nd column is BAM file alias(optional)
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
    help="""Index of column with color levels (1-based); 
        NOTE: LUAD|red -> LUAD while be labeled in plots and red while be the fill color
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
    "-p",
    "--process",
    type=click.IntRange(min=1, max=cpu_count()),
    default=1,
    help="""
    How many cpu to use \b
    """
)
@click.option(
    "--plot-by",
    type=click.INT,
    default=-1,
    help="""
    Index of column with same plot (1-based)
    """,
    show_default=True
)
@click.option(
    "--sep-by-color",
    is_flag=True,
    type=click.BOOL,
    help="whether to plot colors in different plot",
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
    "--title",
    type=click.STRING,
    default=None,
    help="Title",
    show_default=True
)
@click.option(
    "--distance-ratio",
    type=click.FLOAT,
    default=0.3,
    help="distance between transcript label and transcript line",
    show_default=True
)
def line(
        bam,
        event,
        gtf,
        output,
        config,
        indicator_lines,
        share_y,
        no_gene,
        color_factor,
        dpi,
        log,
        process,
        plot_by,
        sep_by_color,
        remove_empty_gene,
        title,
        distance_ratio
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
    :param indicator_lines: draw vertical lines in sashimi to show the spliced sites
    :param share_y: make different plots use same y axis boundary
    :param no_gene: do not show gene id
    :param color_factor: 1-based index, only work with bam list
    :param dpi: output file resolution
    :param log: whether to perform y axis log transform
    :param process:
    :param plot_by:
    :param sep_by_color:
    :param remove_empty_gene:
    :param title
    :param distance_ratio
    :return:
    """
    try:
        log = int(log)
    except ValueError:
        pass

    out_dir = os.path.dirname(os.path.abspath(output))

    try:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    except IOError as err:
        print(err)
        print("Create output directory failed, please check %s" % out_dir)
        exit(err)

    sashimi_plot_settings = parse_settings(config)

    colors = sashimi_plot_settings["colors"]

    # Read bam files
    bam_list, shared_y = prepare_bam_list(bam, color_factor, colors, plot_by=plot_by)

    splice_region = get_sites_from_splice_id(event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        region=splice_region.copy()
    )

    if remove_empty_gene:
        splice_region.remove_empty_transcripts()

    # reads_depth is dict of {BAM: ReadDepth}
    reads_depth = read_reads_depth_from_bam(
        bam_list=bam_list,
        splice_region=splice_region.copy(),
        log=log,
        n_jobs=process
    )

    if share_y:
        assign_max_y(shared_y.values(), reads_depth)

    # shared_y is dict of {group: [BAM...]}
    # then format reads_depth into {group: {BAM: ReadDepth}, group: {BAM: ReadDepth}}
    data = {}
    for key, value in shared_y.items():
        if sep_by_color:
            for val in value:
                if val in reads_depth.keys():
                    k = "{}&%&{}".format(key, val.label)
                    temp = data.get(k, {})
                    temp[val] = reads_depth[val]
                    data[k] = temp
        else:
            data[key] = {x: reads_depth[x] for x in value if x in reads_depth.keys()}

    draw_line_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=data,
        splice_region=splice_region.copy(),
        no_bam=False,
        show_gene=not no_gene,
        dpi=dpi,
        log=log,
        title=title,
        distance_ratio=distance_ratio
    )
