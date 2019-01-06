#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu at 2018.12.19

Inspired by SplicePlot -> mRNAObjects
"""
import os
import re

import pysam
from tqdm import tqdm

from src.data_types import SpliceRegion, ReadDepth, bam_info, GenomicLoci, clean_table_filename
from src.logger import logger


def _is_gtf_(infile):
    u"""
    check if input file is gtf
    :param infile: path to input file
    :return:
    """
    with open(infile) as r:
        for line in r:
            if line.startswith("#"):
                continue

            lines = re.split(r"\s+", line)

            if len(lines) < 8:
                return False

            return bool(
                re.search(
                    r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+",
                    " ".join(lines[8:])
                )
            )


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
            os.remove(infile + ".bai")
            create = True
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


def index_gtf(input_gtf, sort_gtf=False, retry=0):
    u"""
    Created by ygidtu

    Extract only exon tags and keep it clean

    :param input_gtf: path to input gtf file
    :param sort_gtf: Boolean value, whether to sort gtf file first
    :param retry: only try to sort gtf once
    :return path to compressed and indexed bgzipped gtf file
    """
    if not _is_gtf_(input_gtf):
        raise ValueError("gtf file required, %s seems not a valid gtf file" % input_gtf)

    index = False
    output_gtf = input_gtf + ".gz"
    if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
        index = True

    elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf) or \
            os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
        index = True

    # 2018.12.21 used to handle gtf not sorted error
    if sort_gtf and retry > 1:
        raise OSError("Create index for %s failed, and trying to sort it failed too" % input_gtf)
    elif sort_gtf:
        data = []

        logger.info("Sorting %s" % input_gtf)

        old_input_gtf = input_gtf
        input_gtf = re.sub("\.gtf$", "", input_gtf) + ".sorted.gtf"

        output_gtf = input_gtf + ".gz"

        if os.path.exists(input_gtf) and os.path.exists(output_gtf):
            return output_gtf

        try:
            w = open(input_gtf, "w+")
        except IOError as err:
            w = open("/tmp/sorted.gtf")

        with open(old_input_gtf) as r:
            for line in tqdm(r):
                if line.startswith("#"):
                    w.write(line)
                    continue

                lines = line.split()

                data.append(
                    GenomicLoci(
                        chromosome=lines[0],
                        start=lines[3],
                        end=lines[4],
                        strand=lines[6],
                        gtf_line=line
                    )
                )

        for i in sorted(data):
            w.write(i.gtf_line)

        w.close()

    if index:
        try:
            pysam.tabix_index(
                input_gtf,
                preset="gff",
                force=True,
                keep_original=True
            )
        except OSError as err:
            logger.info(err)
            logger.info("Guess gtf needs to be sorted")
            return index_gtf(input_gtf=input_gtf, sort_gtf=True, retry=retry + 1)

    return output_gtf


def read_transcripts(gtf_file, region, retry=0):
    u"""
    Read transcripts from tabix indexed gtf files

    The original function check if the junctions corresponding to any exists exons, I disable this here

    :param gtf_file: path to bgzip gtf files (with tabix index), only ordered exons in this gtf file
    :param region: splice region
    :param retry: if the gtf chromosome and input chromosome does not match. eg: chr9:1-100:+ <-> 9:1-100:+
    :return: SpliceRegion
    """
    if not os.path.exists(gtf_file):
        raise FileNotFoundError("%s not found" % gtf_file)

    splice_region = SpliceRegion(
        chromosome=region.chromosome,
        start=region.start,
        end=region.end,
        strand=region.strand
    )
    try:
        with pysam.Tabixfile(gtf_file, 'r') as gtf_tabix:

            relevant_exons_iterator = gtf_tabix.fetch(
                region.chromosome,
                region.start - 1,
                region.end + 1,
                parser=pysam.asGTF()
            )

            # min_exon_start, max_exon_end, exons_list = float("inf"), float("-inf"),  []
            for line in relevant_exons_iterator:
                splice_region.add_gtf(line)
    except ValueError as err:
        logger.warn(err)

        # handle the mismatch of chromosomes here
        if retry < 2:
            if not region.chromosome.startswith("chr"):
                logger.info("Guess need 'chr'")
                region.chromosome = "chr" + region.chromosome
            else:
                logger.info("Guess 'chr' is redundant")
                region.chromosome = region.chromosome.replace("chr", "")

            return read_transcripts(gtf_file=gtf_file, region=region, retry=retry + 1)
    return splice_region


def read_reads_depth_from_bam(bam_list, splice_region, threshold=0):
    u"""
    read reads coverage info from all bams
    :param bam_list: namedtuple (alias, title, path, label)
    :param splice_region: SpliceRegion
    :param threshold: filter low abundance junctions
    :return: dict {alias, ReadDepth}
    """

    assert isinstance(splice_region, SpliceRegion), "splice_region should be SplcieRegion, not %s" % type(splice_region)

    res = {}
    for bam in bam_list:
        tmp = ReadDepth.determine_depth(
            bam_file_path=bam.path,
            chrm=splice_region.chromosome,
            start_coord=splice_region.start,
            end_coord=splice_region.end,
            threshold=threshold
        )

        tmp.shrink(
            new_low=splice_region.start,
            new_high=splice_region.end
        )

        res[bam] = tmp
    return res


def read_reads_depth_from_count_table(
        count_table,
        splice_region,
        required,
        threshold=0
):
    u"""
    Read junction counts from count_table
    :param count_table:
    :param splice_region:
    :param required:
    :param threshold:
    :return: {label: ReadDepth}
    """

    data = {}
    header = {}
    with open(count_table) as r:
        for line in r:
            lines = line.split()
            # print(lines)
            if not header:
                for i, j in enumerate(lines):
                    header[i] = clean_table_filename(j)
            else:

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
        key = bam_info(path=None, alias=key, label=None, title="")

        res[key] = ReadDepth.create_depth(value, splice_region)

        res[key].shrink(splice_region.start, splice_region.end)

    return res


if __name__ == '__main__':
    pass
