import gzip
import os
import re

import filetype
import pysam
from matplotlib.colors import is_color_like
from src.BamInfo import BamInfo
from src.logger import logger
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
        indicator_lines = {x: 'k' for x in sites}
    elif indicator_lines:
        indicator_lines = [int(i) for i in indicator_lines.split(",")]

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


def prepare_bam_list(bam, color_factor, colors, share_y_by=-1, plot_by=None, barcodes=None):
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

    # load barcode groups
    barcodes_group = {}
    if barcodes:
        with open(barcodes) as r:
            for line in r:
                lines = re.split(r"\t| {2,}", line.strip())

                if len(lines) > 1:
                    key = f"{lines[0]} ({lines[2]})" if len(lines) == 3 else lines[0]
                    temp = barcodes_group.get(lines[0], {})

                    if key not in temp.keys():
                        temp[key] = [lines[1]]
                    else:
                        temp[key].append(lines[1])

                    barcodes_group[lines[0]] = temp

    # load bam list
    shared_y = {}    # {sample group: [BamInfo...]}
    tmp_color = {}
    color_index = 0
    bam_list = []
    with open(bam) as r:
        for line in r:

            lines = re.split(r"\t| {2,}", line.strip())

            if not os.path.exists(lines[0]) and not os.path.isfile(lines[0]):
                logger.warn("wrong input path or input list sep by blank, it should be '\\t'")
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

            temp_barcodes = {}
            if len(lines) > 1:
                if not lines[1] in barcodes_group.keys():
                    temp_barcodes[lines[1]] = None
            else:
                if not is_bam(bam):
                    raise ValueError("%s seem not ba a valid BAM file" % bam)
                
                temp_barcodes = barcodes_group.get(clean_star_filename(bam), None)
                if not temp_barcodes:
                    temp_barcodes[clean_star_filename(bam)] = None

            for alias, barcode in temp_barcodes.items():
                tmp = BamInfo(
                    path=lines[0],
                    alias=alias,
                    title="",
                    label=None,
                    color=tmp_color[color_label.split("|")[0]],
                    barcodes=set(barcode) if barcode else None
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
