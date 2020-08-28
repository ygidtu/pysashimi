#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.01.04

This script contains all the basic data types used by this suite of scripts

For better organization
"""


class GenomicLoci(object):
    u"""
    Created by ygidtu at 2018.12.19

    A base class to handle the position relationships
    """

    __slots__ = [
        "chromosome",
        "start",
        "end",
        "strand",
        "gtf_line"
    ]

    def __init__(self, chromosome, start, end, strand, gtf_line=None):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        """

        self.chromosome = chromosome

        self.start = int(start)
        self.end = int(end)
        self.gtf_line = gtf_line

        if self.end < self.start:
            raise ValueError("End site should bigger than start site, not %d -> %d" % (self.start, self.end))
        if strand == ".":
            strand = "*"
        if strand not in ("+", "-", "*"):
            raise ValueError("strand should be + or -, not %s" % strand)

        self.strand = strand

    def __str__(self):
        u"""
        convert this to string
        :return:
        """
        return "{chromosome}:{start}-{end}:{strand}".format(
            **{
                "chromosome": self.chromosome,
                "start": self.start,
                "end": self.end,
                "strand": self.strand
            }
        )

    def __gt__(self, other):
        u"""
        whether this loci is downstream of other

        Note:
            make sure the wider range is upstream of narrower

            due to the sort of gtf file, therefore the transcript will be ahead of exons
        :param other:
        :return:
        """
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        if self.start != other.start:
            return self.start > other.start

        return self.end < other.end

    def __lt__(self, other):
        u"""
        whether this loci is upstream of other

        Note:
            make sure the wider range is downstream of narrower

            due to the sort of gtf file, therefore the transcript will be ahead of exons
        :param other:
        :return:
        """
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        if self.start != other.start:
            return self.start < other.start

        return self.end > other.end

    def __eq__(self, other):
        u"""
        whether two loci is the same
        :param other:
        :return:
        """
        return self.chromosome == other.chromosome and \
            self.start == other.start and \
            self.end == other.end

    def __add__(self, other):
        u"""
        merge two sites into one
        :param other:
        :return:
        """
        return GenomicLoci(
            chromosome=self.chromosome,
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            strand=self.strand
        )

    def __hash__(self):
        u"""
        generate hash
        :return:
        """
        return hash((self.chromosome, self.start, self.end))

    @property
    def length(self):
        u"""
        :return: int, the length of this loci
        """
        return self.end - self.start

    def is_overlap(self, other):
        u"""
        whether two loci have any overlaps
        :param other: another GenomicLoci and it's children class
        :return: Boolean
        """
        return self.chromosome == other.chromosome and \
            self.start <= other.end and \
            self.end >= other.start

    @classmethod
    def create_loci(cls, string):
        u"""
        Create loci from String
        :param string: chr1:1-100:+
        :return:
        """
        temp = string.split(":")

        if len(temp) == 3:
            chromosome, sites, strand = temp
        elif len(temp) == 2:
            chromosome, sites = temp
            strand = "*"
        else:
            raise ValueError("Failed to decode genomic region: %s" % string)

        start, end = sites.split("-")

        return cls(chromosome, start, end, strand)

