#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu at 2018.12.19

Inspired by SplicePlot -> mRNAObjects
"""
import os
import re
import traceback
from collections import OrderedDict
from multiprocessing import Pool
from typing import Dict

import pysam

from ioutils.utils import clean_star_filename, is_gtf
from src.BamInfo import BamInfo
from src.GenomicLoci import GenomicLoci
from src.ReadDepth import ReadDepth
from src.SpliceRegion import SpliceRegion
from src.logger import logger
from rich.progress import track


def index_gtf(input_gtf, sort_gtf=True, retry=0):
    u"""
    Created by ygidtu

    Extract only exon tags and keep it clean

    :param input_gtf: path to input gtf file
    :param sort_gtf: Boolean value, whether to sort gtf file first
    :param retry: only try to sort gtf once
    :return path to compressed and indexed b-gzipped gtf file
    """
    if input_gtf.endswith(".gz") and os.path.exists(input_gtf + ".tbi"):
        return input_gtf

    if input_gtf is None:
        return None
    gtf = is_gtf(input_gtf)

    if gtf % 10 != 1:
        raise ValueError(f"gtf file required, {input_gtf} seems not a valid gtf file")

    index = False
    if gtf // 10 > 0:
        output_gtf = input_gtf
    else:
        output_gtf = input_gtf + ".gz"
    if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
        index = True

    elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf) or \
            os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
        index = True

    # 2018.12.21 used to handle gtf not sorted error
    if sort_gtf and retry > 1:
        raise OSError(f"Create index for {input_gtf} failed, and trying to sort it failed too")
    elif sort_gtf:
        data = []

        logger.info(f"Sorting {input_gtf}")

        old_input_gtf = input_gtf
        input_gtf = re.sub(r"\.gtf$", "", input_gtf) + ".sorted.gtf"

        output_gtf = input_gtf + ".gz"

        if os.path.exists(input_gtf) and os.path.exists(output_gtf):
            return output_gtf

        w = open(input_gtf, "w+")
        with open(old_input_gtf) as r:
            for line in r:
                if line.startswith("#"):
                    w.write(line)
                    continue

                lines = line.split()

                if len(lines) < 1:
                    continue

                data.append(
                    GenomicLoci(
                        chromosome=lines[0],
                        start=lines[3],
                        end=lines[4],
                        strand=lines[6],
                        gtf_line=line
                    )
                )

        for i in sorted(data, key=lambda x: [x.chromosome, x.start]):
            w.write(i.gtf_line)

        w.close()

    if index:
        logger.info(f"Create index for {input_gtf}")
        try:
            pysam.tabix_index(
                input_gtf,
                preset="gff",
                force=True,
                keep_original=True
            )
        except OSError as err:

            if re.search("could not open", str(err)):
                raise err

            logger.error(err)
            logger.error("Guess gtf needs to be sorted")
            return index_gtf(input_gtf=input_gtf, sort_gtf=True, retry=retry + 1)

    return output_gtf


def read_transcripts(
        gtf_file: str, region: SpliceRegion,
        genome: str = None, retry: int = 0,
        show_id: bool = False, strandless: bool = True
):
    u"""
    Read transcripts from tabix indexed gtf files

    The original function check if the junctions corresponding to any exists exons, I disable this here

    :param gtf_file: path to bgzip gtf files (with tabix index), only ordered exons in this gtf file
    :param region: splice region
    :param retry: if the gtf chromosome and input chromosome does not match. eg: chr9:1-100:+ <-> 9:1-100:+
    :param genome: path to genome fasta file
    :param show_id: whether show gene id or gene name
    :param strandless: True -> show both + and - transcripts;
                        False only show + or -
    :return: SpliceRegion
    """
    if gtf_file is None:
        return region

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"{gtf_file} not found")

    try:
        logger.info(f"Reading from {gtf_file}")

        if genome:
            with pysam.FastaFile(genome) as fa:
                region.sequence = fa.fetch(region.chromosome, region.start - 1, region.end + 1)

        with pysam.Tabixfile(gtf_file, 'r') as gtf_tabix:
            relevant_exons_iterator = gtf_tabix.fetch(
                region.chromosome,
                region.start - 1,
                region.end + 1,
                parser=pysam.asGTF()
            )

            # min_exon_start, max_exon_end, exons_list = float("inf"), float("-inf"),  []
            for line in relevant_exons_iterator:
                try:

                    if not strandless and line.strand != region.strand:
                        continue

                    region.add_gtf(line, show_id=show_id)
                except IndexError as err:
                    logger.error(err)
    except ValueError as err:
        logger.warning(err)

        # handle the mismatch of chromosomes here
        if retry < 2:
            if not region.chromosome.startswith("chr"):
                logger.info("Guess need 'chr'")
                region.chromosome = "chr" + region.chromosome
            else:
                logger.info("Guess 'chr' is redundant")
                region.chromosome = region.chromosome.replace("chr", "")

            return read_transcripts(gtf_file=gtf_file, region=region, retry=retry + 1)

    return region


def __read_from_bam__(args):
    splice_region = args["splice_region"]
    bam = args["bam"]
    threshold = args["threshold"]
    threshold_of_reads = args["threshold_of_reads"]
    log = args["log"]
    reads = args["reads"]
    barcode_tag = args["barcode_tag"]
    strandless = args["strandless"]

    if not splice_region:
        return None

    try:
        if bam.type == "atac":
            tmp = ReadDepth.determine_depth_by_fragments(
                bam=bam,
                chrom=splice_region.chromosome,
                start_coord=splice_region.start,
                end_coord=splice_region.end,
                strand=splice_region.strand,
                strandless=strandless,
                log=log,
            )
        elif bam.type == "depth":
            tmp = ReadDepth.determine_depth_by_depths(
                bam=bam,
                chrom=splice_region.chromosome,
                start_coord=splice_region.start,
                end_coord=splice_region.end,
                strand=splice_region.strand,
                strandless=strandless,
                log=log, groups=args.get("groups", [])
            )
        else:
            tmp = ReadDepth.determine_depth(
                bam=bam,
                chrom=splice_region.chromosome,
                start_coord=splice_region.start,
                end_coord=splice_region.end,
                threshold=threshold,
                threshold_of_reads=threshold_of_reads,
                log=log,
                reads1=reads,
                barcode_tag=barcode_tag,
                required_strand=splice_region.strand if not strandless else None,
                stack=args.get("stack", False),
            )

            tmp.sequence = splice_region.sequence

        if tmp is None:
            return None

        tmp.shrink(
            new_low=splice_region.start,
            new_high=splice_region.end
        )

        return {bam: tmp}
    except (OSError, IOError) as e:
        logger.warning(e)
        return None


def read_reads_depth_from_bam(
        bam_list, splice_region,
        threshold=0, threshold_of_reads=0, log=None, n_jobs=1,
        reads=None, barcode_tag="CB",
        strandless: bool = True,
        stack: bool = False,
        groups: str = ""
) -> Dict[BamInfo, ReadDepth]:
    u"""
    read records coverage info from all bam files
    :param groups:
    :param strandless:
    :param threshold_of_reads:
    :param stack:
    :param barcode_tag:
    :param bam_list: namedtuple (alias, title, path, label)
    :param splice_region: SpliceRegion
    :param threshold: filter low abundance junctions
    :param log
    :param n_jobs
    :param reads: whether to filter out R1/R2
    :return: dict {alias, ReadDepth}
    """
    logger.info("Reading from input files")
    assert isinstance(splice_region, SpliceRegion), f"splice_region should be Splice Region, not {type(splice_region)}"

    res = OrderedDict()

    group_dict = {}
    if groups:
        with open(groups) as r:
            for idx, line in enumerate(r):
                line = line.strip()
                if line not in group_dict.keys():
                    group_dict[line] = []

                group_dict[line].append(idx + 2)

    cmds = []
    new_list = []
    for bam in bam_list:
        if group_dict:
            for key, groups in group_dict.items():
                bam1 = bam.copy()
                bam1.alias = f"{bam1.alias}[{key}]" if bam1.alias else key
                cmds.append({
                    "splice_region": splice_region,
                    "bam": bam1,
                    "threshold": threshold,
                    "threshold_of_reads": threshold_of_reads,
                    "log": log,
                    "reads": reads,
                    "barcode_tag": barcode_tag,
                    "strandless": strandless,
                    "stack": stack,
                    "groups": groups
                })
                new_list.append(bam1)
        else:
            cmds.append({
                "splice_region": splice_region,
                "bam": bam,
                "threshold": threshold,
                "threshold_of_reads": threshold_of_reads,
                "log": log,
                "reads": reads,
                "barcode_tag": barcode_tag,
                "strandless": strandless,
                "stack": stack,
                "groups": groups
            })

    bam_list += new_list

    try:
        with Pool(min(n_jobs, len(bam_list))) as p:
            temp = list(track(p.imap(__read_from_bam__, cmds), total=len(cmds)))

            for x in temp:
                if x is not None:
                    res.update(x)

        # print(len(res))
    except Exception as err:
        logger.error(err)
        traceback.print_exc()
        exit(err)

    if len(res) == 0:
        logger.error("Error reading files, cannot read anything")
        exit(1)

    return {b: res[b] for b in bam_list if b in res.keys()}


def read_reads_depth_from_count_table(
        count_table, splice_region, required,
        colors, threshold=0
) -> Dict[BamInfo, ReadDepth]:
    u"""
    Read junction counts from count_table
    :param count_table: path to count table
    :param splice_region:
    :param required: list of str, which columns are required to draw
    :param threshold: threshold to filter out low abundance junctions
    :param colors: {key: color}
    :return: {label: ReadDepth}
    """

    data = {}
    header = {}
    with open(count_table) as r:
        for line in r:
            lines = line.split()

            if not header:
                for i, j in enumerate(lines):
                    header[i] = clean_star_filename(j)
            else:
                # check file header, to avoid file format error
                if len(header) == len(lines) - 1:
                    logger.info("Change header index due to: Number of headers == number of columns - 1")
                    new_header = {k + 1: v for k, v in header.items()}
                    header = new_header

                for i, j in enumerate(lines):
                    if i == 0:
                        tmp = GenomicLoci.create_loci(lines[0])

                        if not splice_region.is_overlap(tmp):
                            break
                    else:
                        key = header[i]
                        if required:
                            if header[i] in required.keys():
                                key = required[header[i]]
                            else:
                                continue

                        tmp_junctions = data[key] if key in data.keys() else {}

                        if j != "NA" and int(j) >= threshold:
                            tmp_junctions[lines[0]] = int(j)

                        data[key] = tmp_junctions

    res = {}
    for key, value in data.items():

        # customized junctions will introduce string type of key, and list of colors
        # use this try catch to convert key to index to assign colors
        try:
            color = colors[key]
        except TypeError:
            color = colors[len(res) % len(colors)]

        key = BamInfo(
            path=key,
            alias=key,
            label=None,
            title="",
            color=color
        )

        res[key] = ReadDepth.create_depth(value, splice_region)

        res[key].shrink(splice_region.start, splice_region.end)

    return res


if __name__ == '__main__':
    pass
