#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import os
import re
from collections import namedtuple

import click
from openpyxl import load_workbook
from tqdm import tqdm

from src.plot_settings import parse_settings
from src.reading_input import SpliceRegion
from src.reading_input import read_reads_depth_from_bam, read_transcripts, index_gtf, is_bam
from src.sashimi_plot_utils import draw_sashimi_plot

__dir__ = os.path.dirname(os.path.abspath(__file__))
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


def assign_events(events):
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
            res[current] = tmp
            current = i
            tmp = [i]

    if current not in res.keys():
        res[current] = tmp

    return res


def get_merged_event(events, span, indicator_lines):
    u"""
    merged multiple events to huge region, to reduce the IO
    :param events: DataFrame index, output event id
    :param span: see main
    :param indicator_lines: see main
    :return:
    """
    coords = {}
    for e in events:
        tmp = get_sites_from_splice_id(e, span=span, indicator_lines=indicator_lines)
        tmp_key = "%s#%s" % (tmp.chromosome, tmp.strand)

        tmp_list = coords[tmp_key] if tmp_key in coords.keys() else []
        tmp_list.append(tmp)
        coords[tmp_key] = tmp_list

    return coords


def read_info_from_xlsx(xlsx):
    u"""
    Created by ygidtu at 2018.12.25
    read info from input xlsx
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
                tmp[header[i + 1]] = j.value

            data[row[0].value.strip()] = tmp

    bam_list = []
    header = None
    for i in wb.worksheets[1].rows:
        if not header:
            header = True
            continue

        if i[2].value is not None:
            path = i[2].value.strip()

            if not os.path.exists(path):
                raise ValueError("%s not exist" % path)
            elif not is_bam(path):
                raise ValueError("%s is not a BAM" % path)
        else:
            continue

        tmp = bam_info(
            alias=i[0].value if i[0].value is not None else "",
            title=i[1].value if i[1].value is not None else "",
            path=path,
            label=None
        )

        bam_list.append(tmp)

    if not bam_list:
        raise ValueError("No BAM in xlsx")

    return data, bam_list


@click.group(
    context_settings=CONTEXT_SETTINGS,
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
    "--indicator-lines",
    default=None,
    type=click.STRING,
    help="Where to plot additional indicator lines, comma separated int"
)
@click.option(
    "--shared-y",
    default=False,
    if_flag=True,
    type=False,
    help="Whether different sashimi plots shared same y axis"
)
@click.pass_context
@click.version_option(VERSION, message="Current version %(version)s")
def main(
        ctx,
        gtf,
        output,
        config,
        threshold,
        indicator_lines,
        shared_y
):
    u"""
    Welcome

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :param ctx: passed content
    :param bam: list of input BAM files
    :param gtf: path to gtf file
    :param output: path to output file
    :param event: event id, chr:100-200-100-200:+ etc
    :param config: path to config file, default using settings.ini file under this suite of scripts
    :param threshold:
    :param indicator_lines:
    :param shared_y:
    :return:
    """
    ctx.obj['output'] = output
    ctx.obj['config'] = config
    ctx.obj['gtf'] = gtf
    ctx.obj['threshold'] = threshold
    ctx.obj['indicator_lines'] = indicator_lines
    ctx.obj["shared_y"] = shared_y

    pass


@main.command(
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
@click.pass_context
def single(
        ctx,
        bam,
        event,

):
    u"""
    This function is used to plot single sashimi plotting
    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :param ctx: passed parameters from main
    :param bam: list of input BAM files
    :param event: event id, chr:100-200-100-200:+ etc
    :return:
    """
    output = ctx.obj["output"]
    config = ctx.obj["config"]
    gtf = ctx.obj["gtf"]
    indicator_lines = ctx.obj["indicator_lines"]
    threshold = ctx.obj["threshold"]
    shared_y = ctx.obj["shared_y"]

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

    reads_depth = read_reads_depth_from_bam(
        bam_list=bam_list,
        splice_region=splice_region,
        threshold=threshold
    )

    draw_sashimi_plot(
        output_file_path=output,
        settings=sashimi_plot_settings,
        average_depths_dict=reads_depth,
        splice_region=splice_region,
        shared_y=shared_y
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
    "-g",
    "--gtf",
    type=click.Path(exists=True),
    help="Path to gtf file, both transcript and exon tags are necessary"
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
@click.pass_context
def batch(
        ctx,
        input,
        span,
):
    u"""

    This function is used to plot sashimi in batch

    required a specific format of input file

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :param ctx: passed parameters from main
    :param input: input file in specific format
    :param span: str, but must be int or float
    :return:
    """
    output = ctx.obj["output"]
    config = ctx.obj["config"]
    gtf = ctx.obj["gtf"]
    indicator_lines = ctx.obj["indicator_lines"]
    threshold = ctx.obj["threshold"]
    shared_y = ctx.obj["shared_y"]

    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except IOError as err:
        print(err)
        print("Create output directory failed, please check %s" % output)

    data, bam_list = read_info_from_xlsx(input)

    sashimi_plot_settings = parse_settings(config)

    coords = get_merged_event(data.keys(), span=span, indicator_lines=indicator_lines)

    for k, v in tqdm(coords.items()):

        v = assign_events(v)

        for merged, separate in tqdm(v.items()):

            splice_region = read_transcripts(
                gtf_file=index_gtf(input_gtf=gtf),
                region=merged
            )

            reads_depth = read_reads_depth_from_bam(
                bam_list=bam_list,
                splice_region=splice_region,
                threshold=threshold
            )

            for sep in tqdm(separate):
                tmp_reads_depth_dict = {}

                # add label to read_depth
                for i, j in reads_depth.items():
                    tmp_reads_depth_of_bam = j.get_read_depth(sep)

                    try:
                        i = i._replace(label=data[sep.events][i.alias])
                    except KeyError:
                        pass
                    tmp_reads_depth_dict[i] = tmp_reads_depth_of_bam

                draw_sashimi_plot(
                    output_file_path=os.path.join(output, sep.events + ".pdf"),
                    settings=sashimi_plot_settings,
                    average_depths_dict=tmp_reads_depth_dict,
                    splice_region=splice_region.get_region(sep),
                    shared_y=shared_y
                )


if __name__ == '__main__':
    main(obj={})
    pass
