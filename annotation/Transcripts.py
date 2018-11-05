#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
keep coordinate messages of transcripts with its exons
"""


class Transcripts(object):
    u"""
    keep messages of multiples trancripts
    """

    def __init__(self):
        u"""
        init this class
        self.transcripts = {transcript_id: strand}
        self.exons = {transcript_id: [exon coordinate] }
        """
        self.exons = {}
        self.transcripts = {}

    def add_exon(self, transcript, start, end):
        u"""
        add exons to this
        :param transcript: transcript id
        :param start: start site of exon
        :param end: end site of exon
        :return:
        """
        tmp = self.exons[transcript] if transcript in self.exons.keys() else []

        tmp.append(int(start))
        tmp.append(int(end))

        self.exons[transcript] = tmp

    def add_transcript(self, trancript, strand):
        u"""
        add transcript to this
        :param trancript: transcript id
        :param strand: strand
        :return:
        """
        self.transcripts[trancript] = strand

    def __iter__(self):
        u"""
        add iter to this
        :return:
        """
        for transcript, strand in self.transcripts.items():
            yield transcript, strand, self.exons[transcript]

    def __len__(self):
        return len(self.transcripts)
