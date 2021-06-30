#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""

from copy import deepcopy

from src.GenomicLoci import GenomicLoci
from src.Transcript import Transcript


class SpliceRegion(GenomicLoci):
    u"""
    [original description]

    mRNAsObject represents a set of possible mRNA segments

        chrm is the chromosome number
        strand is the strand
        low is the lowest possible coordinate in the set of possible mRNAs
        high is the highest possible coordinate in the set of possible mRNAs
        mRNAs is a list containing possible mRNAs. Each possible mRNA is represented by a list of exons, and
            each exon is represented by a list containing its lowest coordinate and its highest coordinate

    [Now]

    Mainly kept the original functions and structure

    But replace mRNAs[[exon sites], [exon sites]] with SpliceRegion class

    this class is used to collect all the information about the exons and transcripts inside this region
    """

    def __init__(self, chromosome, start, strand, end, sites=None, events=None, ori=None):
        u"""
        init this class
        :param chromosome:  str, chromosome of these
        :param strand: str
        :param start: the very first site of list of exons
        :param end: the very last site of list of exons
        :param sites: list of int, all of the splice sites
        :param events: the source splice events id
        """
        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
        )
        self.sites = self.set_sites(sites)
        self.events = events
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.__transcripts__ = {}  # {transcript_id: namedtuple(gtf proxy of transcript, [gtf proxy of exons])}
        self.ori = ori
        self.sequence = None

        self.__uniq_transcripts__ = set()

    def set_sites(self, sites):
        res = {}

        if sites:
            for s in sites:
                if s not in res.keys():
                    res[s] = "blue"
                else:
                    res[s] = "red"
        return res

    def __str__(self):
        return '{0}:{1}-{2},{3}'.format(
            self.chromosome,
            self.start,
            self.end,
            self.__transcripts__
        )

    def __add__(self, other):
        u"""
        override add of genomic loci
        :param other: SpliceRegion
        :return:
        """
        if self.is_overlap(other):

            return SpliceRegion(
                chromosome=self.chromosome,
                start=min(self.start, other.start),
                end=max(self.end, other.end),
                strand=self.strand,
                sites=self.sites | other.sites
            )

    @property
    def exon_starts(self):
        u"""
        API for extract all exon starts
        :return:
        """
        starts = []
        for i in self.transcripts:
            for j in i.exons:
                starts.append(j.start)
        return sorted(starts)

    @property
    def exon_ends(self):
        u"""
        API for extract all exon ends
        :return:
        """
        ends = []
        for i in self.transcripts:
            for j in i.exons:
                ends.append(j.end)
        return sorted(ends)

    @property
    def exons(self):
        u"""
        API for extract all exons
        :return:
        """
        res = set()

        for i in self.transcripts:
            for j in i.exons:
                res.add(j)

        return sorted(res)

    @property
    def transcripts(self):
        u"""
        convert self.__transcripts__ to list format

        Note: wrapper gtf proxy to list and dict format
            1. there no need to change the code of sashimi plot
            2. shrink the transcript range

        :return: [[{transcript: id, gene: id, exon: []}, {}, {}], [{}]]
        """
        return sorted(
            self.__transcripts__.values(),
            key=lambda x: (x.start, x.end, len(x.exons)),
            reverse=True
        )

    def add_gtf(self, gtf_line, show_id: bool = False):
        u"""
        add new gtf info to this Transcripts class
        :param gtf_line
        :return:
        """
        if isinstance(gtf_line, Transcript):
            gtf_line = deepcopy(gtf_line)

            if gtf_line.start < self.start:
                gtf_line.start = self.start

            if gtf_line.end > self.end:
                gtf_line.end = self.end

            self.__transcripts__[gtf_line.transcript] = gtf_line

        else:
            if gtf_line.feature in ["transcript", "CDS"]:
                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    self.__transcripts__[gtf_line.transcript_id] = Transcript(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start + 1 if gtf_line.start + 1 > self.start else self.start,
                        end=gtf_line.end + 1 if gtf_line.end + 1 < self.end else self.end,
                        strand=gtf_line.strand,
                        transcript_id=gtf_line.transcript_id,
                        gene_id=gtf_line.gene_id,
                        gene = gtf_line.gene_name,
                        transcript = gtf_line.transcript_name,
                        exons=[],
                        show_id = show_id
                    )

            elif gtf_line.feature == "exon":
                if gtf_line.start + 1 >= self.end or gtf_line.end + 1 <= self.start:
                    return

                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    raise ValueError("gtf file not sorted")

                tmp_transcript = self.__transcripts__[gtf_line.transcript_id]

                tmp_transcript.exons.append(
                    GenomicLoci(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start + 1 if gtf_line.start + 1 > tmp_transcript.start else tmp_transcript.start,
                        end=gtf_line.end + 1 if gtf_line.end + 1 < tmp_transcript.end else tmp_transcript.end,
                        strand=gtf_line.strand
                    )
                )

    def get_region(self, genomic):
        u"""
        get smaller splice region
        :return:
        """

        tmp = SpliceRegion(
            chromosome=self.chromosome,
            start=genomic.start,
            end=genomic.end,
            strand=genomic.strand,
            events=genomic.events,
            sites=genomic.sites
        )

        for i in self.transcripts:
            if i.is_overlap(genomic):
                tmp.add_gtf(i)
        return tmp

    def copy(self):
        u"""

        :return:
        """
        temp = SpliceRegion(
            chromosome=self.chromosome,
            start=self.start,
            end=self.end,
            strand=self.strand,
            sites=None,
            events=self.events,
            ori=self.ori
        )

        temp.sites = self.sites

        temp.__transcripts__ = self.__transcripts__
        temp.sequence = self.sequence
        return temp

    def remove_empty_transcripts(self):
        u"""
        Remove transcript without exons
        :return:
        """
        res = {}
        for key, values in self.__transcripts__.items():
            if len(values.exons) > 0:
                res[key] = values

        self.__transcripts__ = res

