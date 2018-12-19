#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Migrated from SplicePlot sashimi_plot_utils
"""

import matplotlib
matplotlib.use('svg')
import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import math
import matplotlib.pyplot as plt

from new_src.transcripts import SpliceRegion


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


def plot_density_single(
        read_depth_object,
        # mRNAs,
        # strand,
        graphcoords,
        graphToGene,
        axvar,
        color='r',
        ymax=None,
        number_junctions=True,
        resolution=.5,
        showXaxis=True,
        nxticks=4,
        font_size=6,
        numbering_font_size=6,
        junction_log_base=10
):
    u"""
    @2018.12.19 remove unnecessary x label
    @2018.12.19 replace junc_comp_function with a new Junction class,
                due to cmp function is removed in Python3

    :param read_depth_object:
    :param graphcoords:
    :param graphToGene:
    :param axvar:
    :param color:
    :param ymax:
    :param number_junctions:
    :param resolution:
    :param showXaxis:
    :param nxticks:
    :param font_size:
    :param numbering_font_size:
    :param junction_log_base:
    :return:
    """
    
    # extract data from read_depth_object
    tx_start = read_depth_object.start

    # @2018.12.19
    # tx_end = read_depth_object.high
    # chrom = read_depth_object.chrm

    wiggle = read_depth_object.wiggle

    jxns = read_depth_object.junctions_dict

    maxheight = max(wiggle)
    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.5 * ymax   

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    prevx = graphcoords[0]
    tmpval = []

    for i in range(len(graphcoords)):
        tmpval.append(wiggle[i])
        if abs(graphcoords[i] - prevx) > resolution:
            compressed_wiggle.append(pylab.mean(tmpval))
            compressed_x.append(prevx)
            prevx = graphcoords[i]
            tmpval = []

    pylab.fill_between(
        compressed_x,
        compressed_wiggle,
        y2=0,
        color=color,
        lw=0
    )

    # @2018.12.19
    # sslists = []
    # for mRNA in mRNAs:
    #     tmp = []
    #     for s, e in mRNA:
    #         tmp.extend([s, e])
    #     sslists.append(tmp)

    # sort the junctions by intron length for better plotting look
    jxns_sorted_list = sorted(jxns.keys())

    current_height = -3 * ymin / 4

    for plotted_count, jxn in enumerate(jxns_sorted_list):
        # leftss, rightss = list(map(int, jxn.split(":")[1].split("-")))

        leftss, rightss = jxn.start, jxn.end

        # @2018.12.19
        # set junctions coordinate here
        # the junction out of boundaries, set the boundaries as coordinate
        ss1 = graphcoords[__get_limited_index__(leftss - tx_start - 1, len(graphcoords))]
        ss2 = graphcoords[__get_limited_index__(rightss - tx_start, len(graphcoords))]

        # ss1, ss2 = graphcoords[leftss - tx_start - 1], graphcoords[

        # mid = (ss1 + ss2) / 2

        # draw junction on bottom
        if plotted_count % 2 == 1:
            pts = [(ss1, 0), (ss1, -current_height), (ss2, -current_height), (ss2, 0)]
            midpt = cubic_bezier(pts, .5)

        # draw junction on top
        else:

            leftdens = wiggle[__get_limited_index__(leftss - tx_start - 1, len(wiggle))]
            rightdens = wiggle[__get_limited_index__(rightss - tx_start, len(wiggle))]

            pts = [(ss1, leftdens),
                   (ss1, leftdens + current_height),
                   (ss2, rightdens + current_height),
                   (ss2, rightdens)]
            midpt = cubic_bezier(pts, .5)

        if number_junctions:
            pylab.text(
                midpt[0],
                midpt[1],
                '{0}'.format(round(jxns[jxn],2)),
                fontsize=numbering_font_size,
                ha='center',
                va='center',
                backgroundcolor='w'
            )

        a = Path(
            pts,
            [
                Path.MOVETO,
                Path.CURVE4,
                Path.CURVE4,
                Path.CURVE4
            ]
        )

        p = PathPatch(
            a,
            ec=color,
            lw=pylab.log(jxns[jxn] + 1) / pylab.log(junction_log_base),
            fc='none'
        )

        axvar.add_patch(p)

    # Format plot
    # ylim(ymin, ymax)
    # axvar.spines['left'].set_bounds(0, ymax)
    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')

    if showXaxis:
        axvar.xaxis.set_ticks_position('bottom')

        # @2018.12.19
        # xlabel('Genomic coordinate (%s), "%s" strand'%(chrom,
        #                                                strand),
        #        fontsize=font_size)

        max_graphcoords = max(graphcoords) - 1

        pylab.xticks(
            pylab.linspace(0, max_graphcoords, nxticks),
            [graphToGene[int(x)] for x in pylab.linspace(0, max_graphcoords, nxticks)],
            fontsize=font_size
        )

    else:
        axvar.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(graphcoords))
    # Return modified axis
    return axvar


# Plot density for a series of bam files.
def plot_density(
        settings,
        event,
        read_depths_dict,
        splice_region,
        # ordered_genotypes_list
):
    u"""
    Several modifications were taken

    1. Due to the change of mRNAObjects, the plot_mRNAs need to be modified
    2. ordered_genotypes_list used to plot the different allele specific situation,
        for now, I changed this to plot multiple BAM files

    :param settings:
    :param event: str, splice event id, only used for subtitle
    :param read_depths_dict:
    :param splice_region:
    :return:
    """

    assert isinstance(splice_region, SpliceRegion)

    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    colors = settings["colors"]
    ymax = settings["ymax"]
    number_junctions = settings["number_junctions"]
    resolution = settings["resolution"]
    junction_log_base = settings["junction_log_base"]
    reverse_minus = settings["reverse_minus"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]
    show_xlabel = settings["show_xlabel"]
    plot_title = settings["plot_title"]
    numbering_font_size = settings["numbering_font_size"]

    # Always show y-axis for read densities for now
    showYaxis = True
    
    # parse mRNA_object to get strand, exon_starts, exon_ends, tx_start, tx_end, chrom
    # strand = mRNA_object.strand
    # chom = mRNA_object.chromosome
    exon_starts = splice_region.exon_starts
    exon_ends = splice_region.exon_ends
    tx_start = splice_region.start
    tx_end = splice_region.end
    mRNAs = splice_region.transcripts
    strand = splice_region.start

    # Get the right scalings
    graphcoords, graphToGene = getScaling(
        tx_start,
        tx_end,
        strand,
        exon_starts,
        exon_ends,
        intron_scale,
        exon_scale,
        reverse_minus
    )

    nfiles = len(list(read_depths_dict.keys()))

    if plot_title is not None and plot_title != '':
        # Use custom title if given
        pylab.suptitle(plot_title, fontsize=10)
    elif plot_title == '':
        pylab.suptitle(event, fontsize=10)
        
    plotted_axes = []

    labels_list = []

    u"""
    Modified by Zhang yiming at 2018.12.19
    
    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """
    # for i, group_genotype in enumerate(ordered_genotypes_list):
    for i, group_genotype in enumerate(read_depths_dict.keys()):
        average_read_depth = read_depths_dict[group_genotype]

        if colors is not None:
            color = colors[i]
        else:
            color = None
        if i < nfiles - 1:
            showXaxis = False 
        else:
            showXaxis = True 

        ax1 = plt.subplot2grid(
            (nfiles + 2, 1),
            (i, 0),
            colspan=1
        )
        
        # Read sample label
        sample_label = group_genotype
        labels_list.append(group_genotype)

        plotted_ax = plot_density_single(
            read_depth_object=average_read_depth,
            # sample_label=sample_label,
            # mRNAs=mRNAs,
            # strand=strand,
            graphcoords=graphcoords,
            graphToGene=graphToGene,
            axvar=ax1,
            # paired_end=False,
            # intron_scale=intron_scale,
            # exon_scale=exon_scale,
            color=color,
            ymax=ymax,
            number_junctions=number_junctions,
            resolution=resolution,
            showXaxis=showXaxis,
            # showYaxis=showYaxis,
            # nyticks=nyticks,
            nxticks=nxticks,
            # show_ylabel=show_ylabel,
            # show_xlabel=show_xlabel,
            font_size=font_size,
            numbering_font_size=numbering_font_size,
            junction_log_base=junction_log_base
        )

        plotted_axes.append(plotted_ax)


    ##
    ## Figure out correct y-axis values
    ##
    ymax_vals = []
    if ymax != None:
        # Use user-given ymax values if provided
        max_used_yval = ymax
    else:
        # Compute best ymax value for all samples: take
        # maximum y across all.
        used_yvals = [curr_ax.get_ylim()[1] for curr_ax in plotted_axes]
        # Round up
        max_used_yval = math.ceil(max(used_yvals))

    # Reset axes based on this.
    # Set fake ymin bound to allow lower junctions to be visible
    fake_ymin = -0.5 * max_used_yval
    universal_yticks = pylab.linspace(
        0,
        max_used_yval,
        nyticks + 1
    )

    # Round up yticks
    universal_ticks = list(map(math.ceil, universal_yticks))
    for sample_num, curr_ax in enumerate(plotted_axes):
        if showYaxis:
            curr_ax.set_ybound(lower=fake_ymin, upper=1.2*max_used_yval)
            curr_yticklabels = []
            for label in universal_yticks:
                if label <= 0:
                    # Exclude label for 0
                    curr_yticklabels.append("")
                else:
                    if label % 1 != 0:
                        curr_yticklabels.append("%.1f" % label)
                    else:
                        curr_yticklabels.append("%d" % label)

            curr_ax.set_yticklabels(
                curr_yticklabels,
                fontsize=font_size
            )

            curr_ax.spines["left"].set_bounds(0, max_used_yval)
            curr_ax.set_yticks(universal_yticks)
            curr_ax.yaxis.set_ticks_position('left')
            curr_ax.spines["right"].set_color('none')

            if show_ylabel:
                y_horz_alignment = 'left'

                curr_ax.set_ylabel(
                    'Depth',
                    fontsize=font_size,
                    va="center",
                    ha=y_horz_alignment,labelpad=10
                )

        else:
            curr_ax.spines["left"].set_color('none')
            curr_ax.spines["right"].set_color('none')
            curr_ax.set_yticks([])
        ##
        ## Plot sample labels
        ##
        sample_color = colors[sample_num]
        # Make sample label y position be halfway between highest
        # and next to highest ytick
        if len(universal_yticks) >= 2:
            halfway_ypos = (universal_yticks[-1] - universal_yticks[-2]) / 2.
            label_ypos = universal_yticks[-2] + halfway_ypos
        else:
            label_ypos = universal_yticks[-1]
        curr_label = labels_list[sample_num]
        curr_ax.text(max(graphcoords), label_ypos,
                     curr_label,
                     fontsize=font_size,
                     va='bottom',
                     ha='right',
                     color=sample_color)
                

    # Draw gene structure
    ax1 = pylab.subplot2grid(
        (nfiles + 2, 1),
        (nfiles, 0),
        colspan=1,
        rowspan=2
    )

    plot_mRNAs(tx_start, mRNAs, graphcoords, reverse_minus)
    pylab.subplots_adjust(hspace=.1, wspace=.7)


def getScaling(
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
    exoncoords = pylab.zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start: exon_ends[i] - tx_start] = 1

    graphToGene = {}
    graphcoords = pylab.zeros((tx_end - tx_start + 1), dtype='f')

    x = 0
    if strand == '+' or not reverse_minus:
        for i in range(tx_end - tx_start + 1):
            graphcoords[i] = x
            graphToGene[int(x)] = i + tx_start
            if exoncoords[i] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    else:
        for i in range(tx_end - tx_start + 1):
            graphcoords[-(i + 1)] = x
            graphToGene[int(x)] = tx_end - i + 1
            if exoncoords[-(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale

    return graphcoords, graphToGene


def plot_mRNAs(
        tx_start,
        mRNAs,
        # strand,
        graphcoords,
        reverse_minus
):
    """
    [original description]
    Draw the gene structure.

    [now]
    Due to I changed the mRNA class, therefore, this function need be modified
    """
    yloc = 0 
    exonwidth = .3
    narrows = 50

    # @2018.12.19
    # the mRNAs is a dict of {transcript: id, gene: id, exon: [Exon, Exon]}
    for mRNA in mRNAs:
        for m in mRNA:
            strand = "+"
            # @2018.12.19
            # s and e is the start and end site of single exon
            for exon in m["exons"]:
                s, e, strand = exon.start, exon.end, exon.strand
                s = s - tx_start
                e = e - tx_start
                x = [
                    graphcoords[s],
                    graphcoords[e],
                    graphcoords[e],
                    graphcoords[s]
                ]
                y = [
                    yloc - exonwidth / 2,
                    yloc - exonwidth / 2,
                    yloc + exonwidth / 2,
                    yloc + exonwidth / 2
                ]
                pylab.fill(x, y, 'k', lw=.5, zorder=20)

            # Draw intron.
            #axhline(yloc, color='k', lw=.5)
            pylab.plot([min(graphcoords), max(graphcoords)], [yloc,yloc], color='k', lw=0.5)

            # Draw intron arrows.
            spread = .2 * max(graphcoords) / narrows
            for i in range(narrows):
                loc = float(i) * max(graphcoords) / narrows
                if strand == '+' or reverse_minus:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
                pylab.plot(x, y, lw=.5, color='k')

            yloc += 1

    pylab.xlim(0, max(graphcoords))
    pylab.ylim(-.5, len(mRNAs) + .5)
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])


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


def draw_sashimi_plot(
        output_file_path,
        settings,
        var_pos,
        average_depths_dict,
        mRNAs_object,
        # ordered_genotypes_list
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

    plt.figure(
        figsize=[settings['width'], settings['height']]
    )

    plot_density(
        settings,                           # plot settings, untouched
        var_pos,                            # splice event id, only used for subtitle
        average_depths_dict,                # reads coverage
        mRNAs_object,                       # Exon and transcript information
        # ordered_genotypes_list              # for now, do not know
    )

    # print("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True
    )
