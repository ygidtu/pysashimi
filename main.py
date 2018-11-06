#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Trying to make a suit of pure python scripts to make sashimi plots
"""
import re
import argparse
import sys
import logging
from annotation.Annotation import Annotation
from annotation.Transcripts import Transcripts
from coverage.Coverage import Coverage
from src.Sashimi import Sashimi
from itertools import zip_longest


class Main(object):
    u"""
    Main program
    """

    def __init__(self):
        u"""
        init this class without any parameter
        """
        self.logger = self.__set_logging__()

        self.parser = argparse.ArgumentParser(description="pysashimi")
        self.__set_parser__()
        self.main()
        pass

    def main(self):

        if len(sys.argv) <= 1:
            self.parser.print_help()
            exit(0)

        try:
            args = self.parser.parse_args(sys.argv[1:])
            print(args)
            tran = Transcripts()
            for i in range(5):
                tran.add_transcript("1", "+")
                tran.add_exon("1", i ** 2 * 10, (i ** 2 + 1) * 10)

            if args.version:
                print("2018.10.28")
                exit(0)

            ref = Annotation(args.reference) if args.reference is not None else None

            if args.which == "index":
                if ref is not None:
                    ref.generate_index()
            elif args.which == "plot":

                chromosome, start, end = re.search(
                    r"(?P<chrom>[\w\.]+):(?P<start>\d+)-(?P<end>\d+)",
                    args.coordinate
                ).groups()

                if args.name is None:
                    names = [None]
                else:
                    names = args.name.split(",")

                coverages = []
                for infile, name in zip_longest(args.infile, names):
                    tmp = Coverage(
                        infile=infile,
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        junctions=args.junctions,
                        name=name
                    )
                    coverages.append(tmp)

                Sashimi(
                    start=start,
                    end=end,
                    transcripts=ref.query(chromosome, start, end) if ref is not None else None,
                    coverages=coverages
                ).save(args.output)

        except argparse.ArgumentError as err:
            self.logger.error(err)
            exit(err)

        pass

    @staticmethod
    def __set_logging__(log_file=None):
        u"""
        set logging style and level
        :return:
        """
        formatter = logging.Formatter('[%(asctime)s] - %(levelname)s - %(message)s')

        # create logger with 'spam_application'
        logger = logging.getLogger('pysashimi')
        logger.setLevel(logging.DEBUG)

        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        # create file handler which logs even debug messages
        if log_file:
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)
        return logger

    def __set_parser__(self):
        self.parser.set_defaults(which="main")
        self.parser.add_argument("-v", "--version", action="store_true", default=False)

        sub_parser = self.parser.add_subparsers(help="Sub-commands")

        generate_index = sub_parser.add_parser("index", help="Generate index for gtf|gff3 file")
        generate_index.set_defaults(which="index")

        generate_index.add_argument(
            "reference", type=str
        )

        # plot subcommand
        plot = sub_parser.add_parser("plot", help="Plot sashimi plot")
        plot.set_defaults(which="plot")

        plot.add_argument(
            "-r",
            "--reference",
            help="Path to gtf|gff3 file",
            type=str,
            # required=True
        )

        plot.add_argument(
            "-c", "--coordinate",
            type=str,
            help="Genomic coordinates, eg: -c 1:100-200",
            required=True
        )

        plot.add_argument(
            "-n", "--name",
            type=str,
            default=None,
            help="The ylabel in final plot ,eg:input1,input2 [default: Names of input BAM/SAM or bigWig files]"
        )

        plot.add_argument(
            "-i", "--infile",
            type=str,
            help="Path to input BAM/SAM/BigWig files, multiple files: -i first -i second",
            required=True,
            action="append"
        )

        plot.add_argument(
            "-s",
            "--strand",
            choices=["NONE", "SENSE", "ANTISENSE", "MATE1_SENSE", "MATE2_SENSE"],
            default="NONE",
            help="Strand specificity [default: %(default)s]"
        )

        plot.add_argument(
            "-j", "--junctions",
            type=str,
            help="Path to bgziped bed files, with tabix index",
            default=None
        )

        plot.add_argument(
            "-t", "--threshold",
            type=int,
            default=0,
            help="Threshold to filter the junctions"
        )

        plot.add_argument(
            "-o", "--output",
            type=str,
            help="Path to output figure",
            required=True
        )


if __name__ == "__main__":
    Main()

