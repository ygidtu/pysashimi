#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
handle gff3 or gtf files
including:
extract transcripts with its exons
generate index
"""
import os
import re
import json
from annotation.Transcripts import Transcripts


class Annotation(object):
    u"""
    I wrote this type of things so many times that I just feel kind of easy?
    """

    def __init__(self, infile):
        u"""
        init this class with input file
        :param infile: path to the input file
        """
        self.__index_gap__ = 5000
        self.infile = os.path.abspath(infile)
        self.__infile_index__ = self.infile + ".idx"

        if not os.path.exists(self.infile) or\
                not os.path.isfile(self.infile):
            raise FileNotFoundError("%s is not existed file" % self.infile)

    def query(self, chromosome, start, end):
        u"""
        read file in
        :return:
        """
        start, end = int(start), int(end)

        file_pointer = 0
        if os.path.exists(self.__infile_index__):
            with open(self.__infile_index__) as r:
                index = json.load(r)
                file_pointer = index[chromosome][start // self.__index_gap__]

        transcript = Transcripts()
        with open(self.infile) as r:
            r.seek(file_pointer)
            for line in r:
                if line.startswith("#"):
                    continue

                lines = line.split("\t")

                if lines[0] < chromosome:
                    continue

                if lines[0] > chromosome:
                    break

                if re.search(r"(transcripts|rna)", lines[2], re.I) and int(lines[3]) > end:
                    break

                if re.search(r"(transcripts|rna)", lines[2], re.I) and int(lines[4]) < start:
                    continue

                if re.search(r"(transcripts|rna)", lines[2], re.I) or \
                        (
                            re.search(r"rna", lines[2], re.I) and
                            "Parent" in lines[8]
                        ):
                    
                    try:
                        transcript.add_transcript(
                            trancript=re.search(
                                    pattern=r"(transcript_id|ID)[= ]\"?(?P<id>\w+)\"?;",
                                    string=lines[8]
                                ).groupdict()["id"],
                            strand=lines[6]
                        )
                    except AttributeError:
                        continue

                if lines[2] == "exon":
                    transcript.add_exon(
                        transcript=re.search(
                                pattern=r"(transcript_id|Parent)[= ]\"?(transcript:)?(?P<id>\w+)\"?;",
                                string=lines[8]
                            ).groupdict()["id"],
                        start=lines[3],
                        end=lines[4]
                    )
        return transcript

    def generate_index(self):
        u"""
        create index for gtf files
        :return:
        """
        # keep index of file pointer, {chrom, [pointer for 0, pointier for 5000, pointer for 10000, ...]}
        data = {}
        with open(self.infile) as r:
            file_pointer = r.tell()
            for line in r:
                if line.startswith("#"):
                    continue

                lines = line.split()

                tmp = data[lines[0]] if lines[0] in data.keys() else []

                index = int(lines[3]) / self.__index_gap__

                if index >= len(tmp):
                    tmp.append(file_pointer)
                    data[lines[0]] = tmp

        with open(self.infile + ".idx", "w+") as w:
            json.dump(data, w, indent=4)

