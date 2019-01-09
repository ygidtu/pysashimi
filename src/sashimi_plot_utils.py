#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Migrated from SplicePlot sashimi_plot_utils

1. change junctions size with fix number; DONE
2. change subtitle to y title; DONE
3. shrink of density not properly working; DONE
4. add transcript id and to mRNAs; DONE
5. tests multiple BAM; DONE
6. display chromosome and strand; DONE
7. Junctions count的过滤筛选; DONE
8. add title, remove title from setting files; DONE
9. add parameter to decide whether use shared y axis
10. fix transcripts display issues
"""
import math

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from src.data_types import SpliceRegion
from src.logger import logger


def __get_limited_index__(num, length):
    u"""
    Created by Zhang yimint at 2018.12.19

    Due to the original author didn't draw any element out of provided range
    So the scripts will through a lot of IndexError

    This function is used to scale that index into the reasonable range

    :param num: current index
    :param length: the list or numpy array length
    :return: int, 0 <= num <= length - 1
    """
    if num < 0:
        return 0

    if num >= length:
        return length - 1

    return num


def cubic_bezier(pts, t):
    """
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = pylab.array(p0)
    p1 = pylab.array(p1)
    p2 = pylab.array(p2)
    p3 = pylab.array(p3)
    return p0 * (1 - t)**3 + 3 * t * p1 * (1 - t) ** 2 + \
        3 * t**2 * (1 - t) * p2 + t**3 * p3


def get_scaling(
        tx_start,
        tx_end,
        strand,
        exon_starts,
        exon_ends,
        intron_scale,
        exon_scale,
        reverse_minus
):
    """
    Compute the scaling factor across various genetic regions.
    """
    exon_coords = pylab.zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exon_coords[exon_starts[i] - tx_start: exon_ends[i] - tx_start] = 1

    graph_to_gene = {}
    graph_coords = pylab.zeros((tx_end - tx_start + 1), dtype='f')

    x = 0
    if strand == '+' or not reverse_minus:
        for i in range(tx_end - tx_start + 1):
            graph_coords[i] = x
            graph_to_gene[int(x)] = i + tx_start
            if exon_coords[i] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    else:
        for i in range(tx_end - tx_start + 1):
            graph_coords[-(i + 1)] = x
            graph_to_gene[int(x)] = tx_end - i + 1
            if exon_coords[-(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale

    return graph_coords, graph_to_gene


def plot_density_single(
        read_depth_object,
        chromosome,
        strand,
        graph_coords,
        graph_to_gene,
        ax_var,
        color='r',
        y_max=None,
        number_junctions=True,
        # resolution=.5,
        show_x_axis=True,
        nx_ticks=4,
        font_size=6,
        numbering_font_size=6,
        no_bam=False

):
    u"""
    @2018.12.19 remove unnecessary x label
    @2018.12.19 replace junc_comp_function with a new Junction class,
                due to cmp function is removed in Python3
    :param chromosome:
    :param strand:
    :param read_depth_object:
    :param graph_coords:
    :param graph_to_gene:
    :param ax_var:
    :param color:
    :param y_max:
    :param number_junctions:
    :param resolution:
    :param show_x_axis:
    :param nx_ticks:
    :param font_size:
    :param numbering_font_size:
    :param no_bam:
    :return:
    """
    
    # extract data from read_depth_object
    tx_start = read_depth_object.start

    # @2018.12.19
    # tx_end = read_depth_object.high
    # chrom = read_depth_object.chrm

    wiggle = read_depth_object.wiggle

    jxns = read_depth_object.junctions_dict

    max_height = read_depth_object.max
    if y_max is None:
        y_max = 1.1 * max_height
    else:
        y_max = y_max
    y_min = -.5 * y_max

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    # prev_x = graph_coords[0]

    u"""
    @2019.01.04
    
    If there is no bam file, use half of y axis as the upper bound of exon 
    
    And draw a white point to maintain the height of the y axis
    
    """
    for i in range(len(graph_coords)):
        compressed_wiggle.append(wiggle[i])
        compressed_x.append(graph_coords[i])

    if no_bam:
        plt.plot(0, max(compressed_wiggle) * 2, color="white")
    else:
        pylab.fill_between(
            compressed_x,
            compressed_wiggle,
            y2=0,
            color=color,
            lw=0
        )

    # sort the junctions by intron length for better plotting look
    jxns_sorted_list = sorted(jxns.keys())

    if not jxns:
        max_junction_count, min_junction_count = 0, 0
    else:
        max_junction_count = max(jxns.values())
        min_junction_count = min(jxns.values())
    junction_count_gap = max_junction_count - min_junction_count

    current_height = -3 * y_min / 4
    for plotted_count, jxn in enumerate(jxns_sorted_list):
        # leftss, rightss = list(map(int, jxn.split(":")[1].split("-")))

        leftss, rightss = jxn.start, jxn.end

        # @2018.12.19
        # set junctions coordinate here
        # the junction out of boundaries, set the boundaries as coordinate
        ss1 = graph_coords[__get_limited_index__(leftss - tx_start - 1, len(graph_coords))]
        ss2 = graph_coords[__get_limited_index__(rightss - tx_start, len(graph_coords))]

        # draw junction on bottom
        if plotted_count % 2 == 1:
            pts = [
                (ss1, 0),
                (ss1, -current_height),
                (ss2, -current_height),
                (ss2, 0)
            ]
            midpt = cubic_bezier(pts, .5)

        # draw junction on top
        else:

            left_dens = wiggle[__get_limited_index__(leftss - tx_start - 1, len(wiggle))]
            right_dens = wiggle[__get_limited_index__(rightss - tx_start, len(wiggle))]

            """
            @2019.01.04
            
            If there is no bam, lower half of y axis as the height of junctions
            """
            if no_bam:
                left_dens /= 2
                right_dens /= 2

            pts = [
                (ss1, left_dens),
                (ss1, left_dens + current_height),
                (ss2, right_dens + current_height),
                (ss2, right_dens)
            ]

            midpt = cubic_bezier(pts, .5)

        if number_junctions:
            t = pylab.text(
                midpt[0],
                midpt[1],
                '{0}'.format(round(jxns[jxn], 2)),
                fontsize=numbering_font_size,
                ha='center',
                va='center',
                backgroundcolor='w'
            )

            # @2018.12.19 transparent background
            t.set_bbox(dict(alpha=0))

        a = Path(
            pts,
            [
                Path.MOVETO,
                Path.CURVE4,
                Path.CURVE4,
                Path.CURVE4
            ]
        )

        """
        @2018.12.26
        scale the junctions line width
        """
        if junction_count_gap > 0:
            line_width = (jxns[jxn] - min_junction_count) / junction_count_gap
        else:
            line_width = 0

        p = PathPatch(
            a,
            ec=color,
            lw=line_width + 0.2,
            fc='none'
        )

        ax_var.add_patch(p)

    # Format plot
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')

    if show_x_axis:
        ax_var.xaxis.set_ticks_position('bottom')

        # @2018.12.19 unnecessary text in figure
        pylab.xlabel(
            'Genomic coordinate (%s), "%s" strand' % (
                chromosome,
                strand
            ),
            fontsize=font_size
        )

        max_graph_coords = max(graph_coords) - 1

        pylab.xticks(
            pylab.linspace(0, max_graph_coords, nx_ticks),
            [graph_to_gene[int(x)] for x in pylab.linspace(0, max_graph_coords, nx_ticks)],
            fontsize=font_size
        )

    else:
        ax_var.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(graph_coords))


def plot_transcripts(
        tx_start,
        transcripts,
        graph_coords,
        reverse_minus,
        font_size,
        show_gene=False,
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
    """
    y_loc = 0
    exon_width = .3

    """
    @2018.12.26
    Maybe I'm too stupid for this, using 30% of total length of x axis as the gap between text with axis
    """
    distance = 0.3 * (max(graph_coords) - min(graph_coords))

    # @2018.12.19
    # @2018.12.21
    # the API of SpliceRegion has changed, the transcripts here should be sorted

    for transcript in transcripts:
        # narrows = math.floor(narrows * (transcript.length / len(graphcoords)))

        # @2018.12.20 add transcript id, based on fixed coordinates
        if show_gene:
            pylab.text(
                x=-1 * distance,
                y=y_loc + 0.15,
                s=transcript.gene,
                fontsize=font_size
            )

            pylab.text(
                x=-1 * distance,
                y=y_loc - 0.25,
                s=transcript.transcript,
                fontsize=font_size
            )
        else:
            pylab.text(
                x=-1 * distance,
                y=y_loc - 0.1,
                s=transcript.transcript,
                fontsize=font_size
            )

        strand = "+"
        # @2018.12.19
        # s and e is the start and end site of single exon
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
            pylab.fill(x, y, 'k', lw=.5, zorder=20)

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
            color='k',
            lw=0.5
        )

        # @2018.12.23 fix intron arrows issues
        # Draw intron arrows.
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

        y_loc += 1

    pylab.xlim(0, max(graph_coords))
    pylab.ylim(-.5, len(transcripts) + .5)
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])


# Plot density for a series of bam files.
def plot_density(
        settings,
        read_depths_dict,
        splice_region,
        show_gene=False,
        title=None,
        share_y=False,
        no_bam=False
):
    u"""
    Several modifications were taken

    1. Due to the change of mRNAObjects, the plot_mRNAs need to be modified
    2. ordered_genotypes_list used to plot the different allele specific situation,
        for now, I changed this to plot multiple BAM files

    :param settings:
    :param read_depths_dict:
    :param splice_region:
    :param show_gene: Boolean value to decide whether show gene id in this graph,
                    used by plot_transcripts()
    :param title: the title of this plot
    :param share_y: whether different sashimi share same y axis
    :param no_bam:
    :return:
    """

    assert isinstance(splice_region, SpliceRegion)

    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    number_junctions = settings["number_junctions"]
    reverse_minus = settings["reverse_minus"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]
    numbering_font_size = settings["numbering_font_size"]

    # parse mRNA_object to get strand, exon_starts, exon_ends, tx_start, tx_end, chrom
    chromosome = splice_region.chromosome
    exon_starts = splice_region.exon_starts
    exon_ends = splice_region.exon_ends
    tx_start = splice_region.start
    tx_end = splice_region.end
    transcripts = splice_region.transcripts
    strand = splice_region.strand

    # Get the right scalings
    graph_coords, graph_to_gene = get_scaling(
        tx_start,
        tx_end,
        strand,
        exon_starts,
        exon_ends,
        intron_scale,
        exon_scale,
        reverse_minus
    )

    n_files = len(read_depths_dict) + ((len(transcripts) // 2) if len(transcripts) > 1 else 1)

    gs = gridspec.GridSpec(n_files, 1)

    if title:
        # Use custom title if given
        pylab.title(title, fontsize=10)

    """
    @ 2018.12.21
    Figure out correct y-axis values

    if share_y is True, compute best ymax value for all samples: take maximum y across all.
    """
    # Round up
    max_used_y_val = math.ceil(max([x.max for x in read_depths_dict.values()]))

    # @2018.12.20 if max_used_y_val is odd, plus one, for better look
    if max_used_y_val % 2 == 1:
        max_used_y_val += 1

    # Reset axes based on this.
    # Set fake ymin bound to allow lower junctions to be visible
    fake_y_min = - 0.5 * max_used_y_val
    universal_y_ticks = pylab.linspace(
        0,
        max_used_y_val,
        nyticks + 1
    )

    """
    @2019.01.07
    calculate the distance between ylabel and y axis
    """
    distance_between_label_axis = max([len(x.alias) for x in read_depths_dict.keys()]) * 2.5

    u"""
    @2018.12.19
    
    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """

    for i, sample_info in enumerate(read_depths_dict.keys()):
        average_read_depth = read_depths_dict[sample_info]

        show_x_axis = (i == len(read_depths_dict) - 1)
        curr_ax = plt.subplot(gs[i, :])

        """
        Re-calculate the y boundary, if do not share same y axis
        @2018.12.20
        if share_y is False, then calculate the best y_limit per axis
        and flush the universal_y_ticks
        """
        if not share_y:
            # Round up
            max_used_y_val = average_read_depth.max
            fake_y_min = - 0.5 * max_used_y_val
            """
            @2018.12.20 if max_used_y_val is odd, plus one, for better look
            @2019.01.07 flush universal_y_ticks if not share_y
            """
            if max_used_y_val % 2 == 1:
                max_used_y_val += 1

            universal_y_ticks = pylab.linspace(
                0,
                max_used_y_val,
                nyticks + 1
            )

        """
        Main body of sashimi
        @2019.01.07 
        """
        # Read sample label
        plot_density_single(
            read_depth_object=average_read_depth,
            chromosome=chromosome,
            strand=strand,
            graph_coords=graph_coords,
            graph_to_gene=graph_to_gene,
            ax_var=curr_ax,
            color=sample_info.color,
            y_max=max_used_y_val,
            number_junctions=number_junctions,
            show_x_axis=show_x_axis,
            nx_ticks=nxticks,
            font_size=font_size,
            numbering_font_size=numbering_font_size,
            no_bam=no_bam
        )

        """
        Indicator lines
        @2018.12.26 
        add indicator lines
        """
        if splice_region.sites:
            for site in splice_region.sites:
                curr_ax.vlines(
                    x=graph_coords[site - tx_start],
                    ymin=0,
                    ymax=max_used_y_val,
                    linestyles="dashed",
                    lw=0.5
                )

        curr_ax.set_ybound(lower=fake_y_min, upper=1.2 * max_used_y_val)
        curr_ax.spines["left"].set_bounds(0, max_used_y_val)
        curr_ax.spines["right"].set_color('none')

        curr_y_tick_labels = []
        for label in universal_y_ticks:
            if label <= 0:
                # Exclude label for 0
                curr_y_tick_labels.append("")
            else:
                if label % 1 != 0:
                    curr_y_tick_labels.append("%.1f" % label)
                else:
                    curr_y_tick_labels.append("%d" % label)

        if not no_bam:

            curr_ax.set_yticklabels(
                curr_y_tick_labels,
                fontsize=font_size
            )

            curr_ax.set_yticks(universal_y_ticks)
            curr_ax.yaxis.set_ticks_position('left')
        else:
            u"""
            @2019.01.04
            
            If there is no bam file, draw a blank y-axis 
            """
            curr_ax.set_yticks([])
            curr_ax.yaxis.set_ticks_position("none")

        """
        Plot y labels
        @2018.12.20 using BAM label as ylabel

        @2019.01.04 change the standards of distance between ylabel and y-axis
        """
        if show_ylabel:
            curr_ax.set_ylabel(
                sample_info.alias,
                fontsize=font_size,
                va="center",
                labelpad=distance_between_label_axis,  # the distance between ylabel with axis
                rotation="horizontal"
            )

        """
        Plot sample labels
        
        @2018.12.20 remove extra text inside sashimi
        @2018.12.25 Add this text back, normally plot title (cell line or tissue) and PSI if exists
        @2018.12.26 use the max_used_yval as y coord
        @2019.01.06 fix bug about sample_num out of colors bound
        """
        if sample_info.label is not None:
            curr_label = sample_info.label
        else:
            curr_label = sample_info.title

        t = curr_ax.text(
            max(graph_coords),
            max_used_y_val,
            curr_label,
            fontsize=font_size,
            va='bottom',
            ha='right',
            color=sample_info.color
        )

        # @218.12.19 set transparent background
        t.set_bbox(dict(alpha=0))

    # Draw gene structure
    """
    @2018.12.26
    add more subplots, based on the number of transcripts
    """
    if len(transcripts) > 0:
        plt.subplot(gs[len(read_depths_dict):, :])

        plot_transcripts(
            tx_start=tx_start,
            transcripts=transcripts,
            graph_coords=graph_coords,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene
        )
    pylab.subplots_adjust(hspace=.1, wspace=.7)


def draw_sashimi_plot(
        output_file_path,
        settings,
        average_depths_dict,
        splice_region,
        share_y,
        no_bam=False,
        show_gene=True,
        dpi=300
):

    """
        draw_sashimi_plot draws the complete sashimi plot

        output_file_path is the file path that the plot will be written to

        settings is a dict containing the settings for the sashimi plot

        var_pos is the location of the SNP, in the format chr1:12345

        average_depths_dict is a dict containing the average read depths by genotype. The keys are the genotypes,
            and the values are ReadDepth objects

        mRNAs_object is an mRNAsObject containing information about the transcript structure

        plot_title is the title of the plot


        return values:
            None. Draws sashimi plot

    """

    assert isinstance(splice_region, SpliceRegion), "splice_region should be SpliceRegion, not %s" % type(splice_region)

    u"""
    @2019.01.04
    If there is no bam, reduce the height of figure
    """
    if no_bam:
        height = settings['height'] * (len(average_depths_dict) + len(splice_region.transcripts)) // 2
    else:
        height = settings['height'] * (len(average_depths_dict) + len(splice_region.transcripts) // 2)

    plt.switch_backend("Agg")

    plt.figure(
        figsize=[
            settings['width'],
            height
        ],
        dpi=dpi
    )

    plot_density(
        settings,                               # plot settings, untouched
        read_depths_dict=average_depths_dict,   # reads coverage
        splice_region=splice_region,            # Exon and transcript information
        show_gene=show_gene,                    # decide whether display gene id in this plot
        share_y=share_y,
        no_bam=no_bam
    )

    logger.info("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True,
        bbox_inches='tight'
    )

