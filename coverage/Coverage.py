#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
extract reads coverage from BAM/SAM or bigWig file
"""
import os
import re
from collections import Counter

import numpy as np
import pyBigWig
import pysam


class Coverage(object):
    u"""
    extract reads coverage from specific region
    """

    def __init__(self, chromosome, start, end, infile, name=None, strand="NONE", junctions=None):
        u"""
        init this
        :param chromosome:
        :param start:
        :param end:
        :param infile: path to BAM/SAM or bigwig file
        :param junctions: if infile is bigWig, the tabix file that record all the junctions needed
        """

        self.infile = os.path.abspath(infile)

        if name is None:
            self.name = os.path.basename(infile)
        else:
            self.name = name

        self.junctions = os.path.abspath(junctions) if junctions else junctions

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)

        try:
            self.coverage, self.introns = self.__read_bam__(strand)
        except IOError:
            self.introns, self.coverage = self.__read_from_bigWig__()

    @staticmethod
    def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):

        # Match
        if CIGAR_op == "M":
            for i in range(pos, pos + CIGAR_len):
                if i < start or i >= end:
                    continue
                ind = i - start
                a[ind] += 1

        # Insertion or Soft-clip
        if CIGAR_op == "I" or CIGAR_op == "S":
            return pos

        # Deletion
        if CIGAR_op == "D":
            pass

        # Junction
        if CIGAR_op == "N":
            don = pos
            acc = pos + CIGAR_len
            if don > start and acc < end:
                junctions[(don, acc)] = junctions.setdefault((don, acc), 0) + 1

        pos = pos + CIGAR_len

        return pos

    @staticmethod
    def flip_read(s, samflag):
        if s == "NONE" or s == "SENSE":
            return 0
        if s == "ANTISENSE":
            return 1
        if s == "MATE1_SENSE":
            if int(samflag) & 64:
                return 0
            if int(samflag) & 128:
                return 1
        if s == "MATE2_SENSE":
            if int(samflag) & 64:
                return 1
            if int(samflag) & 128:
                return 0

    def __read_bam__(self, s):

        # Initialize coverage array and junction dict
        a = {"+": [0] * (self.end - self.start)}
        junctions = {"+": {}}
        if s != "NONE":
            a["-"] = [0] * (self.end - self.start)
            junctions["-"] = {}

        try:
            f = pysam.AlignmentFile(self.infile, require_index=True)
        except IOError:
            pysam.index(self.infile)
            f = pysam.AlignmentFile(self.infile, require_index=True)

        for line in f.fetch(self.chromosome, self.start, self.end):

            samflag, read_start, CIGAR = line.flag, line.reference_start, line.cigarstring

            # Ignore reads with more exotic CIGAR operators
            if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
                continue

            read_strand = [
                "+", "-"][self.flip_read(s, samflag) ^ bool(int(samflag) & 16)]
            if s == "NONE":
                read_strand = "+"

            CIGAR_lens = re.split(r"[MIDNS]", CIGAR)[:-1]
            CIGAR_ops = re.split(r"[0-9]+", CIGAR)[1:]

            pos = read_start

            for n, CIGAR_op in enumerate(CIGAR_ops):
                CIGAR_len = int(CIGAR_lens[n])
                pos = self.count_operator(
                    CIGAR_op,
                    CIGAR_len,
                    pos,
                    self.start,
                    self.end,
                    a[read_strand],
                    junctions[read_strand]
                )
        f.close()

        return list(a.values())[0], list(junctions.values())[0]

    def __read_from_bigWig__(self):
        u"""
        read reads coverage and junctions from bigwig
        :return:
        """
        try:

            bw = pyBigWig.open(self.infile)

            if not bw.isBigWig():
                return None, None

            counts = bw.values(self.chromosome, self.start, self.end)

            if not self.junctions:
                raise FileNotFoundError("junctions file not founded")

            if not os.path.exists(self.junctions + ".tbi"):
                pysam.tabix_index(self.infile)

            f = pysam.TabixFile(self.junctions)
            query = f.fetch(self.chromosome, self.start, self.end)

            introns = []

            for i in query:
                introns.append(
                    [i.query_alignment_start, i.query_alignment_start + i.query_length]
                )

            f.close()

            return Counter(introns), counts
        except RuntimeError:
            return None, None


if __name__ == '__main__':
    test = Coverage("10", 108251945, 108259755, "/Volumes/WD/PacBio/bam/lin_12.bam")
