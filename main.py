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


__dir__ = os.path.dirname(os.path.abspath(__file__))
VERSION = "1.3.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
)
@click.version_option(VERSION, message="Current version %(version)s")
def main():
    u"""
    Welcome

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :return:
    """

    pass


@main.command(
    context_settings=CONTEXT_SETTINGS,
)
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
@click.option(
    "--share-y-by",
    type=click.INT,
    default=-1,
    help="""
    Index of column with share y axis (1-based), Need --share-y\. 
    For example, first 3 bam files use same y axis, and the rest use another
    """,
    show_default=True
)
def normal(
        bam,
        event,
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
        sort_by_color,
        share_y_by
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
    shared_y = {}
    if is_bam(bam):
        bam_list = [
            BamInfo(
                path=bam,
                alias=clean_star_filename(bam),
                title=None,
                label=None,
                color=colors[0]
            )
        ]

    else:
        tmp_color = {}
        color_index = 0
        bam_list = []
        with open(bam) as r:
            for line in r:
                lines = re.split(r"\t| {2,}", line.strip())

                try:
                    color_label = lines[color_factor - 1]
                except IndexError as err:
                    logger.error(err)
                    logger.error("Wrong color factor")
                    logger.error("Your --color-factor is %d" % color_factor)
                    logger.error("Your error line in %s" % lines)

                    exit(err)

                if color_label not in tmp_color.keys():
                    tmp_color[color_label] = colors[color_index % len(colors)]
                    color_index += 1

                if len(lines) > 1:
                    tmp = BamInfo(
                        path=lines[0],
                        alias=lines[1],
                        title="",
                        label=None,
                        color=tmp_color[color_label]
                    )
                else:
                    if not is_bam(bam):
                        raise ValueError("%s seem not ba a valid BAM file" % bam)

                    tmp = BamInfo(
                        path=bam,
                        alias=clean_star_filename(bam),
                        title="",
                        label=None,
                        color=tmp_color[color_label]
                    )
                bam_list.append(tmp)

                if share_y_by > 0:
                    try:
                        tmp = shared_y.get(lines[share_y_by - 1], [])
                        tmp.append(bam_list[-1])
                        shared_y[lines[share_y_by - 1]] = tmp
                    except IndexError as err:
                        logger.error(err)
                        logger.error("Wrong --share-y-by index")
                        logger.error("Your --share-y-by is %d" % share_y_by)
                        logger.error("Your error line in %s" % lines)

                        exit(err)

    if sort_by_color:
        bam_list = sorted(bam_list, key=lambda x: x.color)

    splice_region = get_sites_from_splice_id(event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        region=splice_region.copy()
    )

    reads_depth = read_reads_depth_from_bam(
        bam_list=bam_list,
        splice_region=splice_region.copy(),
        threshold=threshold,
        log=log,
        n_jobs=process
    )

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

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region.copy(),
        no_bam=False,
        show_gene=not no_gene,
        dpi=dpi,
        log=log
    )


@main.command(
    context_settings=CONTEXT_SETTINGS,
)
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


@main.command()
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


if __name__ == '__main__':
    main.add_command(normal)
    main.add_command(no_bam)
    main.add_command(pipeline)
    main()
    pass


