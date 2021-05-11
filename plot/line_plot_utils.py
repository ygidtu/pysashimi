#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.11.20

Make line
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy
from matplotlib import pylab
from src.logger import logger
from src.SpliceRegion import SpliceRegion

from plot.transcripts_plot_utils import plot_transcripts
from plot.utils import get_scaling


def plot_density_single(
        read_depth_object,
        chromosome,
        strand,
        graph_coords,
        ax_var,
        show_x_axis=True,
        nx_ticks=4,
        font_size=6,
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
    :param show_x_axis:
    :param nx_ticks:
    :param font_size:
    :param log:
    :return:
    """

    # Plot all reads depth as line
    for bam, obj in read_depth_object.items():
        wiggle = obj.wiggle

        # Reduce memory footprint by using incremented graphcoords.
        compressed_x = []
        compressed_wiggle = []

        for i in range(len(graph_coords)):
            compressed_wiggle.append(wiggle[i])
            compressed_x.append(graph_coords[i])

        pylab.plot(compressed_x, compressed_wiggle, color=bam.color, lw=1)

    # Format plot
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')

    if show_x_axis:
        ax_var.xaxis.set_ticks_position('bottom')

        # @2018.12.19 unnecessary text in figure

        xlabel = 'Genomic coordinate (%s), "%s" strand' % (
            chromosome,
            strand
        )

        if log in (2, 10):
            xlabel = xlabel + ", y axis is log%d transformed" % log

        pylab.xlabel(
            xlabel,
            fontsize=font_size
        )

        bk = len(graph_coords) // nx_ticks

        linspace, ticks = [], []
        for i in range(0, len(graph_coords), bk):
            linspace.append(graph_coords[i])
            ticks.append(read_depth_object.start + i)

        pylab.xticks(linspace, ticks, fontsize=font_size)

    else:
        ax_var.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(graph_coords))


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
    :return:
    """

    assert isinstance(splice_region, SpliceRegion)

    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    reverse_minus = settings["reverse_minus"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]

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

    """
    @2019.01.07
    calculate the distance between ylabel and y axis
    """
    distance_between_label_axis = max([len(x) for x in read_depths_dict.keys()]) * 2.5

    u"""
    @2018.12.19
    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """

    for i, group in enumerate(read_depths_dict.keys()):
        average_read_depth = read_depths_dict[group]

        show_x_axis = (i == len(read_depths_dict) - 1)
        curr_ax = plt.subplot(gs[i, :])

        if title is not None and i == 0:
            curr_ax.set_title(title, fontsize=10, loc="center")

        # Round up
        max_used_y_val = list(average_read_depth.values())[0].max
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
            show_x_axis=show_x_axis,
            nx_ticks=nxticks,
            font_size=font_size,
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
                group.split("&%&")[0],
                fontsize=font_size,
                va="center",
                labelpad=distance_between_label_axis,  # the distance between ylabel with axis
                rotation="horizontal"
            )

        if "&%&" in group:
            t = curr_ax.text(
                max(graph_coords),
                max_used_y_val,
                group.split("&%&")[1],
                fontsize=font_size,
                va='bottom',
                ha='right',
                color=list(average_read_depth.keys())[0].color
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
            show_gene=show_gene,
            distance_ratio=distance_ratio
        )
    pylab.subplots_adjust(hspace=.15, wspace=.7)


def draw_line_plot(
        output_file_path,
        settings,
        average_depths_dict,
        splice_region,
        no_bam=False,
        show_gene=True,
        dpi=300,
        log=None,
        title=None,
        distance_ratio=0.3
):
    """
        draw_sashimi_plot draws the complete sashimi plot

        output_file_path is the file path that the plot will be written to

        settings is a dict containing the settings for the sashimi plot

        var_pos is the location of the SNP, in the format chr1:12345

        average_depths_dict {group: {BAM: ReadDepth}, group: {BAM: ReadDepth}}

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

    plt.figure(
        figsize=[settings['width'], height],
        dpi=dpi
    )

    plot_density(
        settings,  # plot settings, untouched
        read_depths_dict=average_depths_dict,  # reads coverage
        splice_region=splice_region,  # Exon and transcript information
        show_gene=show_gene,  # decide whether display gene id in this plot
        no_bam=no_bam,
        log=log,
        title=title,
        distance_ratio=distance_ratio
    )

    logger.info("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True,
        bbox_inches='tight'
    )

