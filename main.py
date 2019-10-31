#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import os
import re

import click
from openpyxl import load_workbook

from src.data_types import SpliceRegion, clean_bam_filename, clean_table_filename, bam_info
from src.logger import logger
from src.plot_settings import parse_settings
from src.reading_input import index_gtf
from src.reading_input import is_bam
from src.reading_input import read_reads_depth_from_bam
from src.reading_input import read_reads_depth_from_count_table
from src.reading_input import read_transcripts
from src.sashimi_plot_utils import draw_sashimi_plot

__dir__ = os.path.dirname(os.path.abspath(__file__))
VERSION = "1.2.6"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def get_sites_from_splice_id(string, span=0, indicator_lines=None):
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

    string = string.strip()
    split = string.split("@")

    if not split:
        raise ValueError("Invalid region %s" % string)

    sites = []
    try:
        for i in split:
            if re.search(r"[\w\.]:(\d+-?){2,}:[+-]", i):
                chromosome, tmp_sites, strand = i.split(":")
            elif re.search(r"[\w\.]:(\d+-?){2,}[+-]", i):
                chromosome, tmp_sites = i.split(":")
                tmp_sites, strand = tmp_sites[:-1], tmp_sites[-1]
            else:
                chromosome, tmp_sites = i.split(":")
                strand = "*"

            try:
                for x in tmp_sites.split("-"):
                    sites.append(int(x))
            except ValueError as err:
                logger.error(err)
                logger.error("Contains illegal characters in %s" % string)
                exit(err)
    except ValueError as err:
        logger.error("Invalid format of input region %s" % string)
        exit(err)

    sites = sorted(sites)
    start, end = sites[0], sites[-1]

    try:
        span = int(span)
        start, end = sites[0] - span, sites[-1] + span
    except ValueError:
        try:
            span = float(span) * (sites[-1] - sites[0])
            start, end = sites[0] - span, sites[-1] + span
        except ValueError as err:
            logger.error("Invalid format of span, %s" % str(span))
            exit(err)

    if indicator_lines is True:
        indicator_lines = sites
    elif indicator_lines is not None:
        indicator_lines = [int(x) for x in indicator_lines.split(",")]

    return SpliceRegion(
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand,
        events=string,
        sites=indicator_lines,
        ori=str(string)
    )


def get_merged_event(events, span, indicator_lines):
    u"""
    merged multiple events to huge region, to reduce the IO
    :param events: DataFrame index, output event id
    :param span: see main
    :param indicator_lines: see main
    :return:
    """
    # if they came from same chr and same strand, the put together
    coords = {}
    for e in events:
        tmp = get_sites_from_splice_id(e, span=span, indicator_lines=indicator_lines)
        tmp_key = "%s#%s" % (tmp.chromosome, tmp.strand)

        tmp_list = coords[tmp_key] if tmp_key in coords.keys() else []
        tmp_list.append(tmp)
        coords[tmp_key] = tmp_list

    return coords


def read_info_from_xlsx(xlsx, color_factor, colors):
    u"""
    Created by ygidtu at 2018.12.25
    read info from input xlsx
    :param xlsx: path to input xlsx file
    :param color_factor: 1-based index, to assign colors to different BAM file
    :param colors: list of colors
    :return:
        - data: {splice region: {header: value}}
        - bam_list: list of bam_info
    """
    wb = load_workbook(filename=xlsx)

    data = {}
    header = None
    for row in wb.worksheets[0].rows:
        if not header:
            header = {i: j.value for i, j in enumerate(row)}

            if not header:
                header = True
        else:
            if not row[0].value:
                continue
            tmp = {}
            for i, j in enumerate(row[1:]):
                tmp[str(header[i + 1]).strip()] = j.value

            data[str(row[0].value).strip()] = tmp

    tmp_color = {}
    color_index = 0
    bam_list = []
    header = None
    for i in wb.worksheets[1].rows:
        if not header:
            header = True
            continue

        if i[2].value is not None:
            path = i[2].value.strip()

            # if not os.path.exists(path):
            #     raise ValueError("%s not exist" % path)
            # elif not is_bam(path):
            #     raise ValueError("%s is not a BAM" % path)
        else:
            continue

        try:
            color_label = i[color_factor - 1].value
        except IndexError as err:
            logger.error(err)
            logger.error("Wrong color factor")
            exit(err)

        if color_label not in tmp_color.keys():
            tmp_color[color_label] = colors[color_index % len(colors)]
            color_index += 1

        tmp = bam_info(
            alias=str(i[1].value) if i[1].value is not None else "",
            title=str(str(i[0].value).strip()) if i[0].value is not None else "",
            path=path,
            label=None,
            color=tmp_color[color_label]
        )

        bam_list.append(tmp)

    if not bam_list:
        raise ValueError("No BAM in xlsx")

    return data, bam_list


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
    type=click.Choice(["0", "2", "10"]),
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
        customized_junction
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
    :return:
    """
    log = int(log)
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

    if is_bam(bam):
        bam_list = [
            bam_info(
                path=bam,
                alias=clean_bam_filename(bam),
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
                lines = line.strip().split("\t")

                try:
                    color_label = lines[color_factor - 1]
                except IndexError as err:
                    logger.error(err)
                    logger.error("Wrong color factor")
                    logger.error("Your --color-factor is ", color_factor)
                    logger.error("Your error line in %s" % bam, lines)

                    exit(err)

                if color_label not in tmp_color.keys():
                    tmp_color[color_label] = colors[color_index % len(colors)]
                    color_index += 1

                if len(lines) > 1:
                    tmp = bam_info(
                        path=lines[0],
                        alias=lines[1],
                        title="",
                        label=None,
                        color=tmp_color[color_label]
                    )
                else:
                    if not is_bam(bam):
                        raise ValueError("%s seem not ba a valid BAM file" % bam)

                    tmp = bam_info(
                        path=bam,
                        alias=clean_bam_filename(bam),
                        title="",
                        label=None,
                        color=tmp_color[color_label]
                    )
                bam_list.append(tmp)

    splice_region = get_sites_from_splice_id(event, indicator_lines=indicator_lines)

    splice_region = read_transcripts(
        gtf_file=index_gtf(input_gtf=gtf),
        region=splice_region.copy()
    )

    reads_depth = read_reads_depth_from_bam(
        bam_list=bam_list,
        splice_region=splice_region.copy(),
        threshold=threshold,
        log=log
    )

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
        share_y=share_y,
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
    type=click.Choice(["0", "2", "10"]),
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
        customized_junction
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
    :return:
    """

    log = int(log)

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
                log=log
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
                share_y=share_y,
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
                lines = line.strip().split("\t")

                if not lines:
                    continue

                if len(lines) > 1:
                    required_cols[lines[0]] = lines[1]
                else:
                    required_cols[lines[0]] = clean_table_filename(lines[0])

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
        share_y=share_y,
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


