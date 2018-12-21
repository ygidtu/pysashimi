#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Craeted by Zhang yiming at 2018.12.16

Migrated from SplicePlot GTFline(class)

Extract the exon from GTF file and create index through tabix
"""
import re
import os
import pysam


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


def index_gtf(input_gtf):
    u"""
    Created by Zhang yiming

    Extract only exon tags and keep it clean

    :param input_gtf: path to input gtf file
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

    if index:
        pysam.tabix_index(
            input_gtf,
            preset="gff",
            force=True,
            keep_original=True
        )

    return output_gtf


if __name__ == '__main__':
    pass

