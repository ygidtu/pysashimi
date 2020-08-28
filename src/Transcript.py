#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""

from src.GenomicLoci import GenomicLoci


class Transcript(GenomicLoci):
    u"""
    Created by ygidtu at 2018.12.21

    A class inherit from GenomicLoci, to collect transcript information
    """

    __slots__ = [
        "transcript",
        "gene",
        "exons",
        "is_reads"
    ]

    def __init__(self, chromosome, start, end, strand, transcript_id, gene_id, exons, is_reads=False):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        :param gene_id: str
        :param exons: list of pysam.GTFProxy
        :param is_reads: is flag used by transcript  plot draw
        """

        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )
        self.transcript = transcript_id
        self.gene = gene_id
        self.exons = exons
        self.is_reads = is_reads

    def __str__(self):
        u"""

        :return:
        """

        exons_str = []

        for i in self.exons:
            exons_str.append("{}-{}".format(i.start, i.end))

        return "{}:{}-{}:{} {} {}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.transcript,
            "|".join(exons_str)
        )

    def __hash__(self):

        exons = sorted([str(x.__hash__()) for x in self.exons])

        return hash((self.chromosome, self.start, self.end, self.strand, " ".join(exons)))

