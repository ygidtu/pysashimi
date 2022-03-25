#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.12.25

Merge three gtf to annotation gtf with replace the original gene and mrna id
"""

import re

from collections import OrderedDict
from multiprocessing import Pool

from rich.progress import track
from fire import Fire


important = ["gene_id", "gene_name", "transcript_id", "transcript_name", "exon_id", "exon_name"]


def __extract_information__(line: str) -> OrderedDict:
    u"""
    extract information from gff of gtf file
    :param line: string after column 8 of gff|gtf file, eg: ID=xxx;Name=xxx or gene_id "xxx"; gene_name "xxx"
    :return: dict
    """
    data = OrderedDict({})
    for i in line.split(";"):
        i = i.strip()
        if not i:
            continue

        tmp = re.split(r"[= ]", i)

        tmp_info = tmp[1].strip() if ":" not in tmp[1] else tmp[1].split(":")[1].strip()
        data[tmp[0]] = tmp_info.replace("\"", "")

    return data


def concat_dict_to_string(data: dict) -> str:
    u"""
    将字典中的键值对构成字符串
    """
    res = []

    for i in important:
        if i in data.keys():
            res.append(f"{i} \"{data[i]}\"")

    for k, v in data.items():
        if k not in important:
            res.append(f"{k} \"{v}\"")
    return "; ".join(res)


def read_gtf(path: str, label: str = None) -> list:
    u"""
    read gtf and replace the gene id and transcript id
    :param path:
    :param label:
    :return:
    """
    data = []

    with open(path) as r:
        for line in track(r):
            if line.startswith("#"):
                continue

            if label:
                line = re.sub(r"MSTRG", label, line)

            lines = line.split()

            gene_id = -1
            ref_gene_id = -1
            for i in range(len(lines)):
                if lines[i] == "gene_id":
                    gene_id = i + 1

                elif lines[i] == "ref_gene_id":
                    ref_gene_id = i + 1

            if gene_id > 0 and ref_gene_id > 0:
                lines[gene_id], lines[ref_gene_id] = lines[ref_gene_id], lines[gene_id]

            data.append("{}\t{}".format("\t".join(lines[:8]), " ".join(lines[8:])))

    return data


def main(reference, output, label=None):

    reference = read_gtf(reference, label)
    # reference += read_gtf(w8, "8W")
    # reference += read_gtf(p7, "P7")
    # reference += read_gtf(e14, "E14")

    with open(output, "w+") as w:
        for lines in track(reference):
            w.write("{}\n".format(lines))


if __name__ == '__main__':
    Fire(main)
