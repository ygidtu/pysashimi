#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.01.04

This script contains all the basic data types used by this suite of scripts

For better organization
"""
import os
import re
from collections import namedtuple
from copy import deepcopy

import numpy
import pysam

from src.logger import logger

clean_bam_filename = lambda x: re.sub("([_.]?Aligned.sortedByCoord.out)?.bam", "", os.path.basename(x))
clean_table_filename = lambda x: re.sub("[_.]?SJ.out.tab", "", os.path.basename(x))

bam_info = namedtuple("bam_info", ["alias", "title", "label", "path", "color"])
ax_label = namedtuple("NamedAx", ["Ax", "Label"])   # @2018.12.20 using this to handle the ylabel of different ax


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
        chromosome, sites, strand = string.split(":")

        start, end = sites.split("-")

        return cls(chromosome, start, end, strand)

    # def to_splice_region(self):
    #     u"""
    #     Convert this genomic loci to splice region
    #     :return:
    #     """
    #     return SpliceRegion(
    #         chromosome=self.chromosome,
    #         start=self.start,
    #         end=self.end,
    #         strand=self.strand
    #     )


class Transcript(GenomicLoci):
    u"""
    Created by ygidtu at 2018.12.21

    A class inherit from GenomicLoci, to collect transcript information
    """

    __slots__ = [
        "transcript",
        "gene",
        "exons",
    ]

    def __init__(self, chromosome, start, end, strand, transcript_id, gene_id, exons):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        :param gene_id: str
        :param exons: list of pysam.GTFProxy
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

    def __init__(self, chromosome, start, strand, end, sites=None, events=None):
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
        self.sites = set(sites) if sites else set()
        self.events = events
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.__transcripts__ = {}  # {transcript_id: namedtuple(gtf proxy of transcript, [gtf proxy of exons])}

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

    def add_gtf(self, gtf_line):
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

            if gtf_line.feature == "transcript":
                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    self.__transcripts__[gtf_line.transcript_id] = Transcript(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start if gtf_line.start > self.start else self.start,
                        end=gtf_line.end if gtf_line.end < self.end else self.end,
                        strand=gtf_line.strand,
                        transcript_id=gtf_line.transcript_id,
                        gene_id=gtf_line.gene_id,
                        exons=[]
                    )

            elif gtf_line.feature == "exon":
                if gtf_line.start >= self.end or gtf_line.end <= self.start:
                    return

                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    raise ValueError("gtf file not sorted")

                tmp_transcript = self.__transcripts__[gtf_line.transcript_id]

                tmp_transcript.exons.append(
                    GenomicLoci(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start if gtf_line.start > tmp_transcript.start else tmp_transcript.start,
                        end=gtf_line.end if gtf_line.end < tmp_transcript.end else tmp_transcript.end,
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


class Junction(object):
    u"""
    Created by ygidtu at 2018.12.19

    This is used to collect information of single junction
    And provide relative position comparision
    """

    __slots__ = [
        "chromosome",
        "start",
        "end"
    ]

    def __init__(self, chromosome, start, end):
        u"""
        init this class
        :param chromosome:
        :param start:
        :param end:
        """
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)

        if self.end <= self.start:
            raise ValueError("End site(%d) should bigger than start site(%d)" % (start, end))

    @property
    def length(self):
        u"""
        :return: int, the length of this junction
        """
        return self.end - self.start

    @classmethod
    def create_junction(cls, string):
        u"""
        create Junction from chr1:1-100:+
        :param string: str, chr1:1-100:+ format or chr1:1-100 also work
        :return:
        """
        string = string.split(":")

        chromosome = string[0]
        start, end = string[1].split("-")

        return cls(chromosome=chromosome, start=start, end=end)

    def __hash__(self):
        u"""
        generate hash
        :return:
        """
        return hash((self.chromosome, self.start, self.end))

    def __str__(self):
        u"""
        convert junctions to string
        :return:
        """
        return "{chrom}:{start}-{end}".format(
            **{
                "chrom": self.chromosome,
                "start": self.start,
                "end": self.end
            }
        )

    def __gt__(self, other):
        u"""
        great than
        compare two junction by length
        :param other:
        :return:
        """
        return self.length > other.length

    def __lt__(self, other):
        u"""
        less than
        compare two junction by length
        :param other:a
        :return:
        """
        return self.length < other.length

    def __eq__(self, other):
        u"""
        same length
        :param other:
        :return:
        """
        return self.length == other.length

    def is_overlap(self, other):
        u"""
        whether any overlap with another Junction or GenomicLoci
        :param other:
        :return:
        """

        if self.chromosome != other.chromosome:
            return False

        return self.start < other.end and self.end > other.start

    def is_upstream(self, other):
        u"""
        whether this junction is upstream of other
        :param other:
        :return:
        """
        assert isinstance(other, Junction), "Input should be Junction class"

        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        return self.end < self.start

    def is_downstream(self, other):
        u"""
        whether this junction is downstream of other
        :param other:
        :return:
        """
        assert isinstance(other, Junction), "Input should be Junction class"

        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        return self.start > other.end


class ReadDepth(GenomicLoci):
    u"""
    Migrated from SplicePlot ReadDepth class

    add a parent class to handle all the position comparision
    """

    def __init__(self, chromosome, start, end, wiggle, junctions_dict):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param wiggle: do not know what it is
        :param junctions_dict:
        """
        super().__init__(chromosome, start, end, "+")

        if wiggle is not None:
            assert chromosome is None or self.length + 1 == len(wiggle), "Wiggle length don't correspond to input range"

        self.wiggle = wiggle
        self.junctions_dict = junctions_dict
        self.max = max(self.wiggle)

    @classmethod
    def determine_depth(
            cls,
            bam_file_path,
            chrm,
            start_coord,
            end_coord,
            threshold
    ):
        """
            determine_depth determines the coverage at each base between start_coord and end_coord, inclusive.

            bam_file_path is the path to the bam file used to \
            determine the depth and junctions on chromosome between start_coord and end_coord

        return values:
            depth_vector,
            which is a Numpy array which contains the coverage at each base position between start_coord and end_coord

            spanned_junctions, which is a dictionary containing the junctions supported by reads.
            The keys in spanned_junctions are the
                names of the junctions, with the format chromosome:lowerBasePosition-higherBasePosition
        """
        try:
            with pysam.AlignmentFile(bam_file_path, 'rb') as bam_file:
                relevant_reads = bam_file.fetch(reference=chrm, start=start_coord, end=end_coord)

                depth_vector = numpy.zeros(end_coord - start_coord + 1, dtype='f')
                spanned_junctions = {}

                for read in relevant_reads:

                    # make sure that the read can be used
                    cigar_string = read.cigartuples

                    # each read must have a cigar string
                    if cigar_string is None:
                        continue

                    # read cannot have insertions or deletions
                    contains_indel = False
                    # spans_more_than_one_junction = False
                    for cigar_event in cigar_string:
                        if cigar_event[0] == 1 or cigar_event[0] == 2:
                            contains_indel = True
                            break

                    if contains_indel:
                        continue

                    for index, base_position in enumerate(read.get_reference_positions()):
                        if start_coord <= base_position <= end_coord:
                            depth_vector[base_position - start_coord] += 1

                        # junction spanning case
                        if (index + 1) < len(read.positions) and \
                                base_position + 1 != read.get_reference_positions()[index + 1]:

                            try:
                                junction_name = Junction(
                                    chrm,
                                    base_position + 1,
                                    read.get_reference_positions()[index + 1] + 1
                                )

                                if junction_name not in spanned_junctions:
                                    spanned_junctions[junction_name] = 0

                                spanned_junctions[junction_name] = spanned_junctions[junction_name] + 1
                            except ValueError:
                                continue

            filtered_junctions = {}
            for k, v in spanned_junctions.items():
                if v >= threshold:
                    filtered_junctions[k] = v

            return cls(
                chrm,
                start_coord,
                end_coord,
                depth_vector,
                filtered_junctions
            )
        except IOError:
            logger.error('There is no .bam file at {0}'.format(bam_file_path))
            raise Exception
        except ValueError as err:
            logger.error(bam_file_path)
            logger.error(err)
            exit(err)

    @classmethod
    def create_depth(cls, data, splice_region, depth=10):
        u"""
        Create ReadDepth base on junction dict
        :param data: {junction in string: int}
        :param splice_region:
        :param depth:
        :return:
        """
        junctions_dict = {}
        for key, value in data.items():
            junctions_dict[Junction.create_junction(key)] = value

        depth_vector = numpy.zeros(
            splice_region.end - splice_region.start + 1,
            dtype='f'
        )

        for i in splice_region.exons:
            for j in range(i.start - splice_region.start, i.end - splice_region.start + 1):
                depth_vector[j] = depth

        return cls(
            chromosome=splice_region.chromosome,
            start=splice_region.start,
            end=splice_region.end,
            wiggle=depth_vector,
            junctions_dict=junctions_dict
        )

    @classmethod
    def create_blank(cls):

        """
            create_blank creates an instance of ReadDepth where all of the attributes are None
        """
        return cls("1", 0, 1, None, None)

    def is_invalid(self):
        """
            is_invalid determines whether any of the attributes are None
        """
        return self.junctions_dict is None or (self.chromosome == "1" and self.start == 0 and self.end == 1)

    def shrink(self, new_low, new_high):

        """
            shrink changes the boundaries of the ReadDepth object

            new_low is the new lower genomic coordinate boundary for the ReadDepth object
            new_high is the new upper genomic coordinate boundary for the ReadDepth object

            This method also changes self.wiggle and self.
            junctions_dict so that they only contain data between new_low and new_high

            return value:
                Nothing. Method changes the ReadDepth object
        """
        if new_low < self.start or new_high > self.end:
            raise Exception('New boundaries are not valid, old: %d-%d, new: %d-%d' % (
                self.start,
                self.end,
                new_low,
                new_high
            ))

        # filter through junctions_dict to remove junctions which are no longer in the region
        new_junctions_dict = {}
        for key, value in self.junctions_dict.items():
            ss_low, ss_high = key.start, key.end

            if ss_low >= new_low and ss_high <= new_high:
                new_junctions_dict[key] = value

        self.junctions_dict = new_junctions_dict

        # replace the wiggle
        bottom_index = new_low - self.start
        top_index = bottom_index + (new_high - new_low)

        self.wiggle = self.wiggle[bottom_index:top_index + 1]

        # change the lower and upper bound (last step)
        self.start = new_low
        self.end = new_high

    def __add__(self, other):

        """
            __add__ allows two ReadDepth objects to be added together using the + symbol

            Both self and other must have the same low and high attributes

            return value:
                A new ReadDepth object containing the sum of the two original ReadDepth objects
        """

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chromosome == other.chromosome, 'Cannot add depths from different chromosomes'
        assert self.start == other.start and self.end == other.end, 'Cannot add depths with different start and end'
        new_wiggle = self.wiggle + other.wiggle

        new_junctions_dict = {}

        for key, value in self.junctions_dict.items():
            if key in other.junctions_dict:
                new_junctions_dict[key] = value + other.junctions_dict[key]
            else:
                new_junctions_dict[key] = value

        for key, value in list(other.junctions_dict.items()):
            if key not in self.junctions_dict:
                new_junctions_dict[key] = value

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            new_wiggle,
            new_junctions_dict
        )

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(
            self.chromosome,
            self.start,
            self.end,
            self.wiggle,
            self.junctions_dict
        )

    def divide_by_constant(self, constant):

        """
        divide_by_constant divides self.wiggle and self.junctions_dict by a constant value

        constant is a number

        return value:
            A new ReadDepth object containing the divided values. Method leaves the original ReadDepth object unchanged
        """
        new_wiggle = self.wiggle / constant

        new_junctions_dict = {}
        for key, value in list(self.junctions_dict.items()):
            new_junctions_dict[key] = value * 1.0 / constant

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            new_wiggle,
            new_junctions_dict
        )

    def filter_junctions_dict_for_event(self, splice_event_name):
        """
            filter_junctions_dict_for_event removes all entries frm junctions_dict that cannot possibly be
                involved in the alternative splicing event splice_event_name

            splice_event_name is the name of the alternative splicing event, in the format
                format chr1:17055-17915,chr1:17055-17606,chr1:17055-17233,
                where the numbers represent the genomic coordinates of the splice sites

            return values:
                A new ReadDepth object containing only the relevant junctions
        """
        junction_names_list = splice_event_name.split(',')
        new_junctions_dict = {}
        for junction_name in junction_names_list:
            if junction_name in self.junctions_dict:
                new_junctions_dict[junction_name] = self.junctions_dict[junction_name]

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            self.wiggle,
            new_junctions_dict
        )

    def get_read_depth(self, genomic):
        u"""
        Created by ygidtu at 2018.12.25
        Extract part of data from this class
        :param genomic: a genomic class
        :return:
        """
        if genomic.start == self.start and genomic.end == self.end:
            return self

        targets = [x for x in range(genomic.start - self.start, genomic.end - self.start + 1)]
        wiggle = self.wiggle.take(targets)

        junctions_dict = {}
        for k, v in self.junctions_dict.items():
            if k.is_overlap(genomic):
                junctions_dict[k] = v

        return ReadDepth(
            chromosome=self.chromosome,
            start=genomic.start,
            end=genomic.end,
            wiggle=wiggle,
            junctions_dict=junctions_dict
        )


if __name__ == '__main__':
    pass
