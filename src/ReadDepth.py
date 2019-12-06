#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""

import numpy
import pysam

from scipy.stats import zscore

from conf.logger import logger
from src.GenomicLoci import GenomicLoci
from src.Junction import Junction


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
        threshold,
        log,
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

        :param log:
        """
        try:
            with pysam.AlignmentFile(bam_file_path, 'rb') as bam_file:

                try:
                    relevant_reads = bam_file.fetch(reference=chrm, start=start_coord, end=end_coord)
                except ValueError:
                    if chrm.startswith("chr"):
                        chrm = chrm.replace("chr", "")
                    else:
                        chrm = "chr{}".format(chrm)
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
                    # spans_more_than_one_junction = False
                    for cigar_event in cigar_string:
                        if cigar_event[0] == 1 or cigar_event[0] == 2:
                            continue

                    for index, base_position in enumerate(read.get_reference_positions()):
                        if start_coord <= base_position + 1 <= end_coord:
                            depth_vector[base_position - start_coord + 1] += 1

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

            if log == 10:
                depth_vector = numpy.log10(depth_vector + 1)
            elif log == 2:
                depth_vector = numpy.log2(depth_vector + 1)
            elif log == "zscore":
                depth_vector = zscore(depth_vector)

            return cls(
                chromosome=chrm,
                start=start_coord,
                end=end_coord,
                wiggle=depth_vector,
                junctions_dict=filtered_junctions
            )
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(bam_file_path))
            raise err
        except ValueError as err:
            logger.error(bam_file_path)
            logger.error(err)
            # exit(err)

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
            for j in range(i.start - splice_region.start + 1, i.end - splice_region.start + 2):
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

    def add_customized_junctions(self, other):
        u"""
        Add customized junctions to plot
        :param other:
        :return:
        """
        new_junctions_dict = {}

        for key, value in self.junctions_dict.items():
            if key in other.junctions_dict:
                new_junctions_dict[key] = value + other.junctions_dict[key]
            else:
                new_junctions_dict[key] = value

        for key, value in list(other.junctions_dict.items()):
            if key not in self.junctions_dict:
                new_junctions_dict[key] = value

        self.junctions_dict = new_junctions_dict


if __name__ == '__main__':
    pass
