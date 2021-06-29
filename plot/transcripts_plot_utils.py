#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import math

from matplotlib import pylab

from plot.utils import set_indicator_lines


def plot_transcripts(
    tx_start,
    transcripts,
    graph_coords,
    reverse_minus,
    font_size,
    ymax=None,
    sites=None,
    show_gene=False,
    distance_ratio=0.3,
    color = None
):
    """
    [original description]
    draw the gene structure.

    [now]
    due to i changed the mrna class, therefore, this function need be modified

    :param tx_start: the very start of this plot
    :param graph_coords: numpy array, convert the coord of genome to the coord in this plot
    :param reverse_minus:
    :param transcripts: list of Transcript
    :param font_size: the font size of transcript label
    :param show_gene: Boolean value to decide whether to show gene id in this plot
    :param distance_ratio: distance between transcript label and transcript line
    :param curr_ax: make transcript plot to specific ax
    """
    y_loc = 0
    exon_width = .3

    """
    @2018.12.26
    Maybe I'm too stupid for this, using 30% of total length of x axis as the gap between text with axis
    """
    distance = distance_ratio * (max(graph_coords) - min(graph_coords))
    
    # @2018.12.19
    # @2018.12.21
    # the API of SpliceRegion has changed, the transcripts here should be sorted

    strand = "+"
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
                        y=y_loc + 0.25,
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
            s = s - tx_start
            e = e - tx_start
            x = [
                graph_coords[s],
                graph_coords[e],
                graph_coords[e],
                graph_coords[s]
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
            graph_coords[transcript.start - tx_start],
            graph_coords[transcript.end - tx_start]
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
            max_ = graph_coords[transcript.end - tx_start]
            min_ = graph_coords[transcript.start - tx_start]
            length = max_ - min_
            narrows = math.ceil(length / max(graph_coords) * 50)

            spread = .2 * length / narrows

            for i in range(narrows):
                loc = float(i) * length / narrows + graph_coords[transcript.start - tx_start]
                if strand == '+' or reverse_minus:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [y_loc - exon_width / 5, y_loc, y_loc + exon_width / 5]
                pylab.plot(x, y, lw=.5, color='k')

        y_loc += 1 # if transcript.transcript else .5
    
    pylab.xlim(0, max(graph_coords))
    pylab.ylim(-.5, len(transcripts) + .5)
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])

    set_indicator_lines(tx_start, None, graph_coords, sites, ymax=ymax)
