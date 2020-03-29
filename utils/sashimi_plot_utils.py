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
import numpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from src.SpliceRegion import SpliceRegion
from conf.logger import logger

from utils.transcripts_plot_utils import plot_transcripts


def __get_limited_index__(num, length):
    u"""
    Created by Zhang yimint at 2018.12.19

    Due to the original author didn't draw any element out of provided range
    So the scripts will through a lot of IndexError

    This function is used to scale that index into the reasonable range

    :param num: current index
    :param length: the list or numpy array length
    :return: (int, bool), 0 <= num <= length - 1, and modified or not
    """
    if num < 0:
        return 0, True

    if num >= length:
        return length - 1, True

    return num, False


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

    graph_coords = pylab.zeros((tx_end - tx_start + 1), dtype='f')

    x = 0
    for i in range(0, tx_end - tx_start + 1):
        graph_coords[i] = x

        if exon_coords[i if strand == '+' or not reverse_minus else -(i + 1)] == 1:
            x += 1. / exon_scale
        else:
            x += 1. / intron_scale

    return graph_coords


def plot_density_single(
        read_depth_object,
        chromosome,
        strand,
        graph_coords,
        ax_var,
        color='r',
        y_max=None,
        number_junctions=True,
        show_x_axis=True,
        nx_ticks=4,
        font_size=6,
        numbering_font_size=6,
        no_bam=False,
        log=None
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
    :param show_x_axis:
    :param nx_ticks:
    :param font_size:
    :param numbering_font_size:
    :param no_bam:
    :param log:
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
            lw=0,
            step="post"
        )

    # sort the junctions by intron length for better plotting look
    jxns_sorted_list = sorted(
        jxns.keys(),
        key=lambda x: x.end - x.start,
        reverse=True
    )

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
        ss1_idx, ss1_modified = __get_limited_index__(leftss - tx_start, len(graph_coords))
        ss2_idx, ss2_modified = __get_limited_index__(rightss - tx_start, len(graph_coords))

        u"""
        @2019.01.14

        add two new variables to make it clear which one is index, which one is genomic site 
        """
        ss1 = graph_coords[ss1_idx]
        ss2 = graph_coords[ss2_idx]

        # draw junction on bottom
        if plotted_count % 2 == 0:

            pts = [
                (ss1, 0 if not ss1_modified else -current_height),
                (ss1, -current_height),
                (ss2, -current_height),
                (ss2, 0 if not ss2_modified else -current_height)
            ]
            midpt = cubic_bezier(pts, .5)

        # draw junction on top
        else:

            left_dens = wiggle[ss1_idx]
            right_dens = wiggle[ss2_idx]

            """
            @2019.01.04
            
            If there is no bam, lower half of y axis as the height of junctions
            """
            if no_bam:
                left_dens /= 2
                right_dens /= 2

            pts = [
                (ss1, left_dens if not ss1_modified else left_dens + current_height),
                (ss1, left_dens + current_height),
                (ss2, right_dens + current_height),
                (ss2, right_dens if not ss2_modified else right_dens + current_height)
            ]

            midpt = cubic_bezier(pts, .5)

        if number_junctions:
            t = pylab.text(
                midpt[0], midpt[1],
                '{0}'.format(round(jxns[jxn], 2)),
                fontsize=numbering_font_size,
                ha='center', va='center',
                backgroundcolor='w'
            )

            # @2018.12.19 transparent background
            t.set_bbox(dict(alpha=0))

        a = Path(
            pts,
            [
                Path.MOVETO, Path.CURVE4,
                Path.CURVE4, Path.CURVE4
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
            a, ec=color,
            lw=line_width + 0.2, fc='none'
        )

        ax_var.add_patch(p)

    # Format plot
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')

    if show_x_axis:
        ax_var.xaxis.set_ticks_position('bottom')

        # @2018.12.19 unnecessary text in figure

        xlabel = 'Genomic coordinate (%s), "%s" strand' % (chromosome, strand)

        if log in (2, 10):
            xlabel = xlabel + ", y axis is log%d transformed" % log

        pylab.xlabel(xlabel, fontsize=font_size)

        bk = 1
        if not read_depth_object.sequence:
            bk = len(graph_coords) // nx_ticks

        linspace, ticks = [], []
        for i in range(0, len(graph_coords), bk):
            linspace.append(graph_coords[i])

            temp_txs = tx_start + i
            if read_depth_object.sequence:
                # print(i - read_depth_object.start, nx_ticks)
                if (i - read_depth_object.start) % nx_ticks == 0:
                    temp_txs = "{}\n{}".format(tx_start + i, read_depth_object.sequence[i])
                else:
                    temp_txs = "\n{}".format(read_depth_object.sequence[i])
            ticks.append(temp_txs)

        pylab.xticks(linspace, ticks, fontsize=font_size)

    else:
        ax_var.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(graph_coords))


# Plot density for a series of bam files.
def plot_density(
        settings,
        read_depths_dict,
        splice_region,
        show_gene=False,
        title=None,
        no_bam=False,
        log=None,
        distance_ratio=0.3
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
    :param no_bam:
    :param log: y ticks log transformation or not, 2 -> log2; 10 -> log10，
                do not use `set_yscale` here, because there are junction under the x axis,
                and there coords do not convert into log by matplotlib, so it will cause a lot troubles
    :param distance_ratio: distance between transcript label and transcript line
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
    graph_coords = get_scaling(
        tx_start, tx_end, strand,
        exon_starts, exon_ends,
        intron_scale, exon_scale,
        reverse_minus
    )

    n_files = len(read_depths_dict) + ((len(transcripts) // 2) if len(transcripts) > 1 else 1)

    gs = gridspec.GridSpec(n_files, 1)

    """
    @2019.01.07
    calculate the distance between ylabel and y axis
    """
    distance_between_label_axis = max([len(x) if isinstance(x, str) else len(x.alias) for x in read_depths_dict.keys()]) * 2.5

    u"""
    @2018.12.19

    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """

    for i, sample_info in enumerate(read_depths_dict.keys()):
        average_read_depth = read_depths_dict[sample_info]

        show_x_axis = (i == len(read_depths_dict) - 1)
        curr_ax = plt.subplot(gs[i, :])

        if title is not None and i == 0:
            curr_ax.set_title(title, fontsize=10, loc="center")

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
            ax_var=curr_ax,
            color=sample_info.color,
            y_max=max_used_y_val,
            number_junctions=number_junctions,
            show_x_axis=show_x_axis,
            nx_ticks=nxticks,
            font_size=font_size,
            numbering_font_size=numbering_font_size,
            no_bam=no_bam,
            log=log
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

        u"""
        The y ticks are formatted here
        @2019.03.31 add little check here to make sure the y axis shows the real value
        """
        curr_y_tick_labels = []

        for label in universal_y_ticks:
            if log in (2, 10):
                label = numpy.power(log, label)

            if label <= 0:
                # Exclude label for 0
                curr_y_tick_labels.append("")
            else:
                curr_y_tick_labels.append("%.1f" % label if label % 1 != 0 else "%d" % label)

                # if log in (2, 10):
                #     curr_y_tick_labels[-1] = r"$\mathregular{{" + str(log) + "}^{" + curr_y_tick_labels[-1] + "}}$"

        if not no_bam:
            curr_ax.set_yticks(universal_y_ticks)
            curr_ax.set_yticklabels(
                curr_y_tick_labels,
                fontsize=font_size
            )

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
            curr_label = ""

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
        plt.subplot(gs[len(read_depths_dict):, :]) # + 1 if splice_region.sequence else len(read_depths_dict)

        plot_transcripts(
            tx_start=tx_start,
            transcripts=transcripts,
            graph_coords=graph_coords,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene,
            distance_ratio=distance_ratio
        )

        if splice_region.sites:
            for site in splice_region.sites:
                plt.vlines(
                    x=graph_coords[site - tx_start],
                    ymin=-0.5,
                    ymax=len(transcripts) - .5,
                    linestyle="dashed",
                    lw=0.5
                )

    pylab.subplots_adjust(hspace=.15, wspace=.7)


def draw_sashimi_plot(
        output_file_path,
        settings,
        average_depths_dict,
        splice_region,
        no_bam=False,
        show_gene=True,
        dpi=300,
        log=None,
        distance_ratio=0.3,
        title=None
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

        if splice_region.sequence:
            height += (settings["height"] * .2)

    plt.figure(figsize=[settings['width'], height], dpi=dpi)

    plot_density(
        settings,                               # plot settings, untouched
        read_depths_dict=average_depths_dict,   # reads coverage
        splice_region=splice_region,            # Exon and transcript information
        show_gene=show_gene,                    # decide whether display gene id in this plot
        no_bam=no_bam,
        log=log,
        distance_ratio=distance_ratio,
        title=title
    )

    logger.info("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True,
        bbox_inches='tight'
    )

