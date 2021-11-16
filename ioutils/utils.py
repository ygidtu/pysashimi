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
    return re.sub("[_.]?(SJ.out.tab|Aligned.sortedByCoord.out)?.bam", "", os.path.basename(x))


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


def load_barcode(path: str) -> dict:
    barcodes = {}
    with open(path) as r:
        for line in r:
            lines = re.split(r"\t| {2,}", line.strip())

            lines[0] = clean_star_filename(lines[0])

            if len(lines) >= 3:
                key = lines[2]
                temp = barcodes.get(lines[0], {})

                if key not in temp.keys():
                    temp[key] = [lines[1]]
                else:
                    temp[key].append(lines[1])

                barcodes[lines[0]] = temp

    return barcodes


def load_colors(bam: str, barcodes: str, color_factor: str, colors):

    res = {}

    if color_factor and re.search("^\\d+$", color_factor):
        color_factor = int(color_factor) - 1

    if color_factor and not isinstance(color_factor, int):
        with open(color_factor) as r:
            for line in r:
                line = line.strip().split("\t")
                res[line[0]] = line[1]

    try:
        with open(bam) as r:
            for idx, line in enumerate(r):
                line = line.strip().split()
                key = line[1] if len(line) > 1 else clean_star_filename(line[0])

                if key not in res.keys():
                    if not isinstance(color_factor, int):
                        res[key] = colors[idx % len(colors)]
                    else:
                        if len(line) <= color_factor:
                            logger.error("--color-factor must <= number of columns from " + bam)
                            exit(1)
                        res[key] = line[color_factor].upper()
                        if "|" in res[key]:
                            res[key] = res[key].split("|")[1]
    except Exception as err:
        logger.error("please check the input file, including: .bai index", err)
        exit(0)

    if barcodes:
        temp = set()
        with open(barcodes) as r:
            for line in r:
                line = line.strip().split()
                temp.add(line[-1])

        for idx, group in enumerate(sorted(temp)):
            if group not in res.keys():
                res[group] = colors[idx % len(colors)]

    return res


def prepare_bam_list(bam, color_factor, colors, share_y_by=-1, plot_by=None, barcodes=None, is_atac = False):
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
    barcodes_group = load_barcode(barcodes) if barcodes else {}

    colors = load_colors(bam, barcodes, color_factor, colors)
     
    # load bam list
    shared_y = {}    # {sample group: [BamInfo...]}
    bam_list = {}
    with open(bam) as r:
        for line in r:

            lines = re.split(r"\t| {2,}", line.strip())

            if not os.path.exists(lines[0]) and not os.path.isfile(lines[0]):
                logger.warn("wrong input path or input list sep by blank, it should be '\\t'")
                continue

            if not is_atac and not is_bam(lines[0]):
                raise ValueError("%s seem not ba a valid BAM file" % lines[0])

            temp_barcodes = barcodes_group.get(
                clean_star_filename(lines[0]), 
                barcodes_group.get(lines[1], {}) if len(lines) > 1 else {}
            )

            if not barcodes_group:
                key = lines[1] if len(lines) > 1 else clean_star_filename(lines[0])
                temp_barcodes[key] = None

            for alias, barcode in temp_barcodes.items():
                tmp = BamInfo(
                    path=lines[0],
                    alias=alias,
                    title="",
                    label=alias,
                    color=colors[alias],
                    barcodes=set(barcode) if barcode else None
                )
       
                if alias not in bam_list.keys():
                    bam_list[alias] = tmp
                else:
                    bam_list[alias] += tmp

                if plot_by is not None:

                    if plot_by < 0:
                        tmp = shared_y.get("", [])
                        tmp.append(alias)
                        shared_y[""] = tmp
                    else:
                        try:
                            tmp = shared_y.get(lines[plot_by - 1], [])
                            tmp.append(alias)
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
                        tmp.append(alias)
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

    return list(bam_list.values()), {i: [bam_list[k] for k in j] for i, j in shared_y.items()}
