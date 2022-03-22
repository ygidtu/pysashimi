#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import math

from matplotlib import pylab

from plot.utils import add_additional_background
from src.SpliceRegion import SpliceRegion


def plot_transcripts(
        transcripts,
        region: SpliceRegion,
        reverse_minus,
        font_size,
        show_gene=False,
        distance_ratio=0.3,
        color=None,
):
    """
    [original description]
    draw the gene structure.

    [now]
    due to i changed the mrna class, therefore, this function need be modified

    :param region:
    :param color:
    :param reverse_minus:
    :param transcripts: list of Transcript
    :param font_size: the font size of transcript label
    :param show_gene: Boolean value to decide whether to show gene id in this plot
    :param distance_ratio: distance between transcript label and transcript line
    """
    y_loc = 0
    exon_width = .3

    """
    @2018.12.26
    Maybe I'm too stupid for this, using 30% of total length of x axis as the gap between text with axis
    """
    distance = distance_ratio * (max(region.graph_coords) - min(region.graph_coords))

    # @2018.12.19
    # @2018.12.21
    # the API of SpliceRegion has changed, the transcripts here should be sorted
    for transcript in transcripts:
        # narrows = math.floor(narrows * (transcript.length / len(graphcoords)))
        strand = transcript.strand
        # @2018.12.20 add transcript id, based on fixed coordinates
        if transcript.transcript:
            if show_gene and transcript.gene:
                if transcript.show_id:
                    pylab.text(
                        x=-1 * distance,
                        y=y_loc + 0.25,
                        s=transcript.gene_id,
                        fontsize=font_size
                    )

                    pylab.text(
                        x=-1 * distance,
                        y=y_loc - 0.25,
                        s=transcript.transcript_id,
                        fontsize=font_size
                    )
                else:
                    pylab.text(
                        x=-1 * distance,
                        y=y_loc,
                        s=transcript.gene + " | " + transcript.transcript,
                        fontsize=font_size
                    )
            else:
                pylab.text(
                    x=-1 * distance,
                    y=y_loc - 0.1,
                    s=transcript.transcript,
                    fontsize=font_size
                )

        # @2018.12.19
        # s and e is the start and end site of single exon
        # print(transcript)
        for exon in transcript.exons:
            s, e, strand = exon.start, exon.end, exon.strand

            x = [
                region.get_relative(s),
                region.get_relative(e),
                region.get_relative(e),
                region.get_relative(s)
            ]
            y = [
                y_loc - exon_width / 2,
                y_loc - exon_width / 2,
                y_loc + exon_width / 2,
                y_loc + exon_width / 2
            ]
            pylab.fill(x, y, 'k' if not color else color, lw=.5, zorder=20)

        # @2018.12.21
        # change the intron range
        # Draw intron.
        intron_sites = [
            region.get_relative(transcript.start),
            region.get_relative(transcript.end),
        ]
        pylab.plot(
            intron_sites,
            [y_loc, y_loc],
            color='k' if not color else color,
            lw=0.5
        )

        # @2018.12.23 fix intron arrows issues
        # Draw intron arrows.
        if not transcript.is_reads:
            length = region.get_relative(transcript.end) - region.get_relative(transcript.start)
            narrows = math.ceil(length / max(region.graph_coords) * 50)

            spread = .2 * length / narrows

            for i in range(narrows):
                loc = float(i) * length / narrows + region.get_relative(transcript.start)
                if strand == '+' or reverse_minus:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [y_loc - exon_width / 5, y_loc, y_loc + exon_width / 5]
                pylab.plot(x, y, lw=.5, color='k')

        y_loc += 1  # if transcript.transcript else .5

    pylab.xlim(0, max(region.graph_coords))
    pylab.ylim(-.5, len(transcripts) + .5)
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])

    if region:
        add_additional_background(region)


if __name__ == '__main__':
    pass
