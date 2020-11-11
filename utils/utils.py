import re
import os
import gzip

import pysam
import filetype

from openpyxl import load_workbook
from matplotlib.colors import is_color_like

from conf.logger import logger
from src.BamInfo import BamInfo
from src.SpliceRegion import SpliceRegion


def clean_star_filename(x):
    u"""
    if input file is STAR SJ.out.tab or STAR bam then remove the SJ.out.tab
    :param x:
    :return:
    """
    return re.sub("[_.]?(SJ.out.tab|Aligned.sortedByCoord.out?.bam)", "", os.path.basename(x))


def is_gtf(infile):
    u"""
    check if input file is gtf
    :param infile: path to input file
    :return:
    """
    if infile is None:
        return False

    is_gtf = 0
    try:
        if filetype.guess_mime(infile) == "application/gzip":
            is_gtf += 10
            r = gzip.open(infile, "rt")
        else:
            r = open(infile)
    except TypeError as err:
        logger.error("failed to open %s", infile)
        exit(err)

    for line in r:
        if line.startswith("#"):
            continue

        lines = re.split(r"\s+", line)

        if len(lines) < 8:
            break

        if re.search(
            r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+",
            " ".join(lines[8:])
        ):
            is_gtf += 1

        break

    r.close()

    return is_gtf


def is_bam(infile):
    u"""
    check if input file is bam or sam file
    :param infile: path to input file
    :return: Boolean
    """

    try:
        create = False
        if not os.path.exists(infile + ".bai"):
            create = True
        elif os.path.getctime(infile + ".bai") < os.path.getctime(infile):
            try:
                os.remove(infile + ".bai")
                create = True
            except PermissionError as err:
                logger.warn(err)
                create = False
        else:
            try:
                with pysam.AlignmentFile(infile) as r:
                    r.check_index()
            except ValueError:
                create = True

        if create:
            logger.info("Creating index for %s" % infile)
            pysam.index(infile)
        return True

    except pysam.utils.SamtoolsError:
        return False


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
            if re.search(r"[\w.]:(\d+-?){2,}:[+-]", i):
                chromosome, tmp_sites, strand = i.split(":")
            elif re.search(r"[\w.]:(\d+-?){2,}[+-]", i):
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

        tmp = BamInfo(
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


def assign_max_y(shared_y, reads_depth, batch = False):
    u"""
    assign max y for input files

    :param shared_y: [[group1], [group2]]
    :param reads_depth: output from read_reads_depth_from_bam
    :return:
    """

    if len(shared_y) == 0:

        max_ = max([x.max for x in reads_depth.values()])

        for v in reads_depth.values():
            v.max = max_
    elif batch:
        max_ = max([reads_depth[x].max for x in shared_y if x in reads_depth.keys()])
        for j in shared_y:
            reads_depth[j].max = max_
    else:

        for i in shared_y:

            max_ = max([reads_depth[x].max for x in i if x in reads_depth.keys()])

            for j in i:
                if j in reads_depth.keys():
                    reads_depth[j].max = max_


def prepare_bam_list(bam, color_factor, colors, share_y_by=-1, plot_by=None):
    u"""
    Prepare bam list
    :return:
    """
    if is_bam(bam):
        return [
            BamInfo(
                path=bam,
                alias=clean_star_filename(bam),
                title=None,
                label=None,
                color=colors[0]
            )
        ], {}

    shared_y = {}    # {sample group: [BamInfo...]}
    tmp_color = {}
    color_index = 0
    bam_list = []
    with open(bam) as r:
        for line in r:

            lines = re.split(r"\t| {2,}", line.strip())

            if not os.path.exists(lines[0]) and not os.path.isfile(lines[0]):
                logger.warning("wrong input path or input list sep by blank, it should be '\\t'")
                continue

            try:
                color_label = lines[color_factor - 1]
            except IndexError as err:
                logger.error(err)
                logger.error("Wrong color factor")
                logger.error("Your --color-factor is %d" % color_factor)
                logger.error("Your error line in %s" % lines)

                exit(err)

            # 如果文件中指定的为特定的颜色，则直接使用该颜色
            if is_color_like(color_label):
                tmp_color[color_label.split("|")[0]] = color_label
            elif "|" in color_label and is_color_like(color_label.split("|")[-1]):
                tmp_color[color_label.split("|")[0]] = color_label.split("|")[-1]
            elif color_label.split("|")[0] not in tmp_color.keys():
                tmp_color[color_label.split("|")[0]] = colors[color_index % len(colors)]
                color_index += 1

            if len(lines) > 1:
                tmp = BamInfo(
                    path=lines[0],
                    alias=lines[1],
                    title="",
                    label=None,
                    color=tmp_color[color_label.split("|")[0]]
                )
            else:
                if not is_bam(bam):
                    raise ValueError("%s seem not ba a valid BAM file" % bam)

                tmp = BamInfo(
                    path=bam,
                    alias=clean_star_filename(bam),
                    title="",
                    label=None,
                    color=tmp_color[color_label.split("|")[0]]
                )
            bam_list.append(tmp)

            if plot_by is not None:
                bam_list[-1].label = color_label.split("|")[0]

                if plot_by < 0:
                    tmp = shared_y.get("", [])
                    tmp.append(bam_list[-1])
                    shared_y[""] = tmp
                else:
                    try:
                        tmp = shared_y.get(lines[plot_by - 1], [])
                        tmp.append(bam_list[-1])
                        shared_y[lines[plot_by - 1]] = tmp
                    except IndexError as err:
                        logger.error(err)
                        logger.error("Wrong --plot-by index")
                        logger.error("Your --plot-by is %d" % plot_by)
                        logger.error("Your error line in %s" % lines)

                        exit(err)

            elif share_y_by > 0:
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

    if len(bam_list) == 0:
        logger.error("Cannot find any input bam file, please check the bam path or the input list")
        exit(1)

    return bam_list, shared_y
