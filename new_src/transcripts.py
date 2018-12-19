#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang yiming at 2018.12.19

Inspired by SplicePlot -> mRNAObjects
"""
import os
import re

import numpy
import pysam


class GenomicLoci(object):
    u"""
    Created by Zhang yiming at 2018.12.19

    A base class to handle the position relationships
    """

    __slots__ = [
        "chromosome",
        "start",
        "end",
        "strand"
    ]

    def __init__(self, chromosome, start, end, strand):
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

        if self.end < self.start:
            raise ValueError("End site should bigger than start site, not %d -> %d" % (self.start, self.end))

        if strand not in ("+", "-"):
            raise ValueError("strand should be + or -, not %s" % strand)

        self.strand = strand

    def __gt__(self, other):
        u"""
        whether this loci is downstream of other
        :param other:
        :return:
        """
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        if self.start != other.start:
            return self.start > other.start

        return self.end > other.end

    def __lt__(self, other):
        u"""
        whether this loci is upstream of other
        :param other:
        :return:
        """
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end

    def __eq__(self, other):
        u"""
        whether two loci is the same
        :param other:
        :return:
        """
        return self.chromosome == other.chromosome and \
               self.start == other.start and \
               self.end == other.end

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


class Exon(GenomicLoci):
    u"""
    Migrated from SplicePlot Exons class
    """

    __slots__ = [
        "transcript",
        "gene"
    ]

    def __init__(self, chromosome, start, end, strand, transcript, gene):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        :param transcript: transcript id, that exon belong to
        :param gene: gene id that exon belong to
        """
        super().__init__(chromosome, start, end, strand)

        self.transcript = transcript
        self.gene = gene

    def __add__(self, other):
        u"""
        merge two Exon into one dict
        :param other: another exon with same transcript with this one
        :return:
        """
        tmp = None
        if self.transcript == other.transcript:
            tmp = self.to_dict()
            tmp["exons"].append(other)
        return tmp

    def __str__(self):
        u"""
        convert this to str
        :return: str
        """
        return "{chrom}:{start}-{end}:{strand}\t{transcript}\t{gene}".format(
            **{
                "chrom": self.chromosome,
                "start": self.start,
                "end": self.end,
                "strand": self.strand,
                "transcript": self.transcript,
                "gene": self.gene
            }
        )

    @staticmethod
    def __extract_info_from_annotation__(gtf_line):
        u"""
        extract chromosome, start, end, strand, transcript id and gene id from gtf line
        :param gtf_line:
        :return:
        """
        data = {}
        for i in gtf_line.split(";"):
            i = i.strip()
            if not i:
                continue

            tmp = re.split(r"[= ]", i)

            tmp_info = tmp[1].strip() if ":" not in tmp[1] else tmp[1].split(":")[1].strip()
            data[tmp[0]] = tmp_info.replace("\"", "")

        return data

    @classmethod
    def create_from_gtf(cls, line):
        u"""
        [original description]

        create_from_gtf creates an Exon object from a single line of a gtf file

        line is a single line of text (a string) from a gtf file

        [now]
        never thought to create a class that way

        :param line: gtf line
        :return:
        """

        info = line.strip('\n').split("\t")
        data = cls.__extract_info_from_annotation__(info[8])

        return cls(
            chromosome=info[0],
            start=int(info[3]),
            end=int(info[4]),
            strand=info[6],
            transcript=data["transcript_id"] if "transcript_id" in data.keys() else "",
            gene=data["gene_id"] if "gene_id" in data.keys() else ""
        )

    def determine_proportion_covered(self, read_depth):
        u"""

        [now]
        still do not understand what this do

        :param read_depth:
        :return:
        """
        assert isinstance(read_depth, ReadDepth), "read_depth should be ReadDepth class, not %s" % type(read_depth)

        if read_depth.chromosome != self.chromosome or \
                self.start < read_depth.start or \
                self.end > read_depth.end:
            return 0

        # determine the top and bottom indices
        bottom_index = self.start - read_depth.start
        top_index = bottom_index + self.length

        bases_covered = 0
        for item in range(bottom_index, top_index + 1):
            if read_depth.wiggle[item] > 0:
                bases_covered += 1.0

        return bases_covered / (top_index + 1 - bottom_index)

    def determine_average_coverage(self, read_depth):
        return self.determine_proportion_covered(read_depth) / (self.end - self.start + 1.0)

    def to_dict(self):
        u"""
        convert Exon class to dict
        :return: dict {transcript: id, gene: id, exons: [Exon, Exon]}
        """
        return {"transcript": self.transcript, "gene": self.gene, "exons": [self]}


class SpliceRegion(object):
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

    def __init__(self, chromosome, start, strand, end):
        u"""
        init this class
        :param chromosome:  str, chromosome of these
        :param strand: str
        :param start: the very first site of list of exons
        :param end: the very last site of list of exons

        """

        self.chromosome = chromosome
        # self.strand = strand
        self.start = start
        self.end = end
        self.strand = strand
        self.__transcripts__ = {}
        # :param transcripts: a list of dict object, [{transcript: id, gene: id, exons: [Exon, Exon]}]

        self.__exon_starts__ = []
        self.__exon_ends__ = []

    def __str__(self):
        return '{0}:{1}-{2},{3}'.format(
            self.chromosome,
            self.start,
            self.end,
            self.__transcripts__
        )

    @property
    def exon_starts(self):
        u"""
        ordered exon starts
        :return:
        """
        return sorted(self.__exon_starts__)

    @property
    def exon_ends(self):
        u"""
        ordered exon ends
        :return:
        """
        return sorted(self.__exon_ends__)

    @property
    def min_start(self):
        u"""
        get the very first start site of exons
        :return: int
        """
        if self.exon_starts:
            return min(self.exon_starts[0], self.start)
        else:
            return self.start

    @property
    def max_end(self):
        u"""
        get the very last end site of exons
        :return: int
        """
        if self.exon_ends:
            return max(self.exon_ends[-1], self.end)
        else:
            return self.end

    @property
    def transcripts(self):
        u"""
        convert self.__transcripts__ to dict format

        :return: [[{transcript: id, gene: id, exon: []}, {}, {}], [{}]]
        """
        data = []
        for i in self.__transcripts__.values():
            data.append([x.to_dict() for x in i])
        return data

    def add_exon(self, exon):
        u"""
        add new exon into this Transcripts class
        :param exon: object of Exon
        :return:
        """
        assert isinstance(exon, Exon), "exon should be an object of Exon, not %s" % type(exon)

        if exon.start >= self.end or exon.end <= self.start:
            return

        tmp = self.__transcripts__.get(exon.transcript)

        if tmp is None:
            tmp = []

        tmp.append(exon)

        self.__transcripts__[exon.transcript] = tmp
        self.__exon_starts__.append(exon.start)
        self.__exon_ends__.append(exon.end)


class Junction(object):
    u"""
    Created by Zhang yiming at 2018.12.19

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
            assert chromosome is None or self.length + 1 == len(wiggle), 'Wiggle, lower bound, and upper bound do not correspond'

        self.wiggle = wiggle
        self.junctions_dict = junctions_dict

    @classmethod
    def determine_depth(cls, bam_file_path, chrm, start_coord, end_coord):
        """
            determine_depth determines the coverage at each base between start_coord and end_coord, inclusive.

            bam_file_path is the path to the bam file used to determine the depth and junctions on chrm between start_coord and end_coord

            return values:
                depth_vector, which is a Numpy array which contains the coverage at each base position between start_coord and end_coord
                spanned_junctions, which is a dictionary containing the junctions supported by reads. The keys in spanned_junctions are the
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

            return cls(chrm, start_coord, end_coord, depth_vector, spanned_junctions)
        except IOError:
            print('There is no .bam file at {0}'.format(bam_file_path))
            raise Exception

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

            This method also changes self.wiggle and self.junctions_dict so that they only contain data between new_low and new_high

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
        assert self.start == other.start and self.end == other.end, 'Cannot add depths with different start and end points'
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


def read_transcripts(gtf_file, chromosome, start, end, strand):
    u"""
    Read transcripts from tabix indexed gtf files

    The original function check if the junctions corresponding to any exists exons, I disable this here

    :param gtf_file: path to bgzip gtf files (with tabix index), only ordered exons in this gtf file
    :param chromosome: chromosome
    :param start: start site
    :param end: end site
    :return: SpliceRegion
    """

    if not os.path.exists(gtf_file):
        raise FileNotFoundError("%s not found" % gtf_file)

    splice_region = SpliceRegion(
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand
    )

    with pysam.Tabixfile(gtf_file, 'r') as gtf_tabix:
        relevant_exons_iterator = gtf_tabix.fetch(
            chromosome,
            start - 1,
            end + 1
        )

        # min_exon_start, max_exon_end, exons_list = float("inf"), float("-inf"),  []
        for line in relevant_exons_iterator:

            exon = Exon.create_from_gtf(line=line)

            if exon.start < start:
                exon.start = start

            if exon.end > end:
                exon.end = end

            splice_region.add_exon(exon)
    return splice_region


def read_reads_depth(bam_list, splice_region, alias=None):
    u"""
    read reads coverage info from all bams
    :param bam_list: list of path to BAM files
    :param splice_region: SpliceRegion
    :param alias: dict {BAM path: BAM alias}
    :return: dict {alias, ReadDepth}
    """

    assert isinstance(splice_region, SpliceRegion), "splice_region should be SplcieRegion, not %s" % type(splice_region)

    res = {}
    for bam in bam_list:
        tmp = ReadDepth.determine_depth(
            bam,
            splice_region.chromosome,
            splice_region.start,
            splice_region.end
        )

        tmp.shrink(
            new_low=splice_region.min_start,
            new_high=splice_region.max_end
        )

        # reduce unnecessary characters
        label = re.sub(r"[_.]SJ.out.tab", "", os.path.basename(bam))
        if alias and bam in alias.keys():
            label = alias[bam]

        res[label] = tmp
    return res




