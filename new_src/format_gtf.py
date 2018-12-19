#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Craeted by Zhang yiming at 2018.12.16

Migrated from SplicePlot GTFline(class)

Extract the exon from GTF file and create index through tabix
"""
import os
import pysam
from tqdm import tqdm


def format_gtf(input_gtf, output_gtf):
    u"""
    Created by Zhang yiming

    Extract only exon tags and keep it clean

    :param input_gtf: path to input gtf file
    :param output_gtf: path to output gtf file
    """

    format_ = False
    if not os.path.exists(output_gtf):
        format_ = True

    else:
        with open(output_gtf) as r:
            for line in tqdm(r):
                if line.startswith("#"):
                    continue

                lines = line.split()
                if lines[2] != "exon":
                    format_ = True
                    break

    if format_:
        data = []
        with open(output_gtf, "w+") as w:
            with open(input_gtf) as r:
                for line in tqdm(r):
                    if line.startswith("#"):
                        continue

                    lines = line.split()

                    if lines[2] == "exon":
                        data.append([lines[0], int(lines[3]), int(lines[4]), line])

                for i in sorted(data, key=lambda x: (x[0], x[1], x[2])):
                    w.write(i[3])

    index = False
    if not os.path.exists(output_gtf + ".gz") and not os.path.exists(output_gtf + ".tbi"):
        index = True

    elif os.path.getctime(output_gtf + ".gz") < os.path.getctime(output_gtf) or \
            os.path.getctime(output_gtf + ".gz") < os.path.getctime(output_gtf):
        index = True

    if index:
        pysam.tabix_index(output_gtf, preset="gff", force=True)


if __name__ == '__main__':
    import sys
    format_gtf(sys.argv[1], sys.argv[2])

