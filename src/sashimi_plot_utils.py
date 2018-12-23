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
from collections import namedtuple
import matplotlib
matplotlib.use('svg')
from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import math
import matplotlib.pyplot as plt

from src.reading_input import SpliceRegion, logger


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


ax_label = namedtuple("NamedAx", ["Ax", "Label"])   # @2018.12.20 using this to handle the ylabel of different ax


def plot_density_single(
        read_depth_object,
        chromosome,
        strand,
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
            compressed_wiggle.append(sum(tmpval) / len(tmpval))
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

            leftdens = wiggle[__get_limited_index__(leftss - tx_start - 1, len(wiggle))]
            rightdens = wiggle[__get_limited_index__(rightss - tx_start, len(wiggle))]

            pts = [
                (ss1, leftdens),
                (ss1, leftdens + current_height),
                (ss2, rightdens + current_height),
                (ss2, rightdens)
            ]

            midpt = cubic_bezier(pts, .5)

        if number_junctions:
            t = pylab.text(
                midpt[0],
                midpt[1],
                '{0}'.format(round(jxns[jxn],2)),
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

        p = PathPatch(
            a,
            ec=color,
            lw=math.log(jxns[jxn] + 1) / math.log(junction_log_base),
            fc='none'
        )

        axvar.add_patch(p)

    # Format plot
    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')

    if showXaxis:
        axvar.xaxis.set_ticks_position('bottom')

        # @2018.12.19 unnecessary text in figure
        pylab.xlabel(
            'Genomic coordinate (%s), "%s" strand'%(
                chromosome,
                strand
            ),
            fontsize=font_size
        )

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
        read_depths_dict,
        splice_region,
        show_gene=False,
        title=None,
        shared_y=False
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
    :param shared_y: whether different sashimi share same y axis
    :return:
    """

    assert isinstance(splice_region, SpliceRegion)

    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    colors = settings["colors"]
    number_junctions = settings["number_junctions"]
    resolution = settings["resolution"]
    junction_log_base = settings["junction_log_base"]
    reverse_minus = settings["reverse_minus"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]
    numbering_font_size = settings["numbering_font_size"]

    # Always show y-axis for read densities for now
    showYaxis = True
    
    # parse mRNA_object to get strand, exon_starts, exon_ends, tx_start, tx_end, chrom
    chromosome = splice_region.chromosome
    exon_starts = splice_region.exon_starts
    exon_ends = splice_region.exon_ends
    tx_start = splice_region.start
    tx_end = splice_region.end
    transcripts = splice_region.transcripts
    strand = splice_region.strand

    # Get the right scalings
    graphcoords, graphToGene = get_scaling(
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

    if title:
        # Use custom title if given
        pylab.title(title, fontsize=10)
        
    plotted_axes = []

    """
    @ 2018.12.21
    Figure out correct y-axis values

    if shared_y is True, compute best ymax value for all samples: take maximum y across all.
    """
    used_yvals = [x.max for x in read_depths_dict.values()]

    # Round up
    max_used_yval = math.ceil(max(used_yvals))

    # @2018.12.20 if max_used_yval is odd, plus one, for better look
    if max_used_yval % 2 == 1:
        max_used_yval += 1

    u"""
    @2018.12.19
    
    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """
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
        plotted_ax = plot_density_single(
            read_depth_object=average_read_depth,
            # sample_label=sample_label,
            chromosome=chromosome,
            strand=strand,
            graphcoords=graphcoords,
            graphToGene=graphToGene,
            axvar=ax1,
            color=color,
            ymax=max_used_yval,
            number_junctions=number_junctions,
            resolution=resolution,
            showXaxis=showXaxis,
            nxticks=nxticks,
            font_size=font_size,
            numbering_font_size=numbering_font_size,
            junction_log_base=junction_log_base
        )

        # @2018.12.16 change ax to [ax, label]
        plotted_axes.append(ax_label(Ax=plotted_ax, Label=group_genotype))

    # Reset axes based on this.
    # Set fake ymin bound to allow lower junctions to be visible
    fake_ymin = -0.5 * max_used_yval
    universal_yticks = pylab.linspace(
        0,
        max_used_yval,
        nyticks + 1
    )

    # Round up yticks
    for sample_num, curr_ax in enumerate(plotted_axes):
        if showYaxis:

            # @2018.12.20
            # if shared_y is False, then calculate the best ylimit per axis
            if not shared_y:
                used_yvals = [curr_ax.Ax.get_ylim()[1] for curr_ax in plotted_axes]

                # Round up
                max_used_yval = math.ceil(max(used_yvals))

                # @2018.12.20 if max_used_yval is odd, plus one, for better look
                if max_used_yval % 2 == 1:
                    max_used_yval += 1

            curr_ax.Ax.set_ybound(lower=fake_ymin, upper=1.2 * max_used_yval)

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

            curr_ax.Ax.set_yticklabels(
                curr_yticklabels,
                fontsize=font_size
            )

            curr_ax.Ax.spines["left"].set_bounds(0, max_used_yval)
            curr_ax.Ax.set_yticks(universal_yticks)
            curr_ax.Ax.yaxis.set_ticks_position('left')
            curr_ax.Ax.spines["right"].set_color('none')

            if show_ylabel:

                # @2018.12.20 using BAM label as ylabel
                curr_ax.Ax.set_ylabel(
                    curr_ax.Label,
                    fontsize=font_size,
                    va="center",
                    labelpad=(len(curr_ax.Label) // 5 + 1) * 10,            # the distance between ylabel with axis
                    rotation="horizontal"
                )

        else:
            curr_ax.Ax.spines["left"].set_color('none')
            curr_ax.Ax.spines["right"].set_color('none')
            curr_ax.Ax.set_yticks([])

        '''
         # @2018.12.20 remove extra text inside sashimi 
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
        t = curr_ax.Ax.text(
            max(graphcoords),
            label_ypos,
            curr_label,
            fontsize=font_size,
            va='bottom',
            ha='right',
            color=sample_color
        )

        # @218.12.19 set transparent background
        t.set_bbox(dict(alpha=0))
        '''

    # Draw gene structure
    ax1 = pylab.subplot2grid(
        (nfiles + 2, 1),
        (nfiles, 0),
        colspan=1,
        rowspan=2
    )

    plot_transcripts(
        tx_start=tx_start,
        transcripts=transcripts,
        graphcoords=graphcoords,
        reverse_minus=reverse_minus,
        font_size=font_size,
        show_gene=show_gene,
        # label_x_axis=graphcoords[tx_start] - 300
    )
    pylab.subplots_adjust(hspace=.1, wspace=.7)


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


def plot_transcripts(
        tx_start,
        transcripts,
        graphcoords,
        reverse_minus,
        font_size,
        show_gene=False,
        # label_x_axis=-300
):
    """
    [original description]
    draw the gene structure.

    [now]
    due to i changed the mrna class, therefore, this function need be modified

    :param tx_start: the very start of this plot
    :param tx_end: the very end of this plot
    :param graphcoords: numpy array, convert the coord of genome to the coord in this plot
    :param reverse_minus:
    :param font_size: the font size of transcript label
    :param show_gene: Boolean value to decide whether to show gene id in this plot
    :param label_x_axis: the position of transcript ids in x axis
    """
    yloc = 0
    exonwidth = .3
    narrows = 50

    # @ 2018.12.21
    # 11258 -> 0.025
    # 21258 -> 0.06
    # this distance is calculated by simple line regression
    distance = 0.005 * math.ceil(len(graphcoords) / 2000)

    # @2018.12.19
    # @2018.12.21
    # the API of SpliceRegion has changed, the transcripts here should be sorted

    for transcript in transcripts:

        # @2018.12.20 add transcript id, based on fixed coordinates
        if show_gene:
            pylab.text(
                x=min(graphcoords) - distance * len(graphcoords),
                y=yloc + 0.1,
                s=transcript.gene,
                fontsize=font_size
            )

            pylab.text(
                x=min(graphcoords) - distance * len(graphcoords),
                y=yloc - 0.2,
                s=transcript.transcript,
                fontsize=font_size
            )
        else:
            pylab.text(
                x=min(graphcoords) - distance * len(graphcoords),
                y=yloc - 0.1,
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

        # @2018.12.21
        # change the intron range
        # Draw intron.
        intron_sites = [graphcoords[transcript.start - tx_start], graphcoords[transcript.end - tx_start]]
        pylab.plot(
            intron_sites,
            [yloc, yloc],
            color='k',
            lw=0.5
        )

        # @2018.12.23 fix intron arrows issues
        # Draw intron arrows.
        max_ = graphcoords[transcript.end - tx_start]

        if max_ < narrows * 3:
            narrows = max_ // 10

        spread = .2 * max_ / narrows

        for i in range(narrows):
            loc = float(i) * max_ / narrows + graphcoords[transcript.start - tx_start]
            if strand == '+' or reverse_minus:
                x = [loc - spread, loc, loc - spread]
            else:
                x = [loc + spread, loc, loc + spread]
            y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
            pylab.plot(x, y, lw=.5, color='k')

        yloc += 1

    pylab.xlim(0, max(graphcoords))
    pylab.ylim(-.5, len(transcripts) + .5)
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
        average_depths_dict,
        splice_region,
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

    assert isinstance(splice_region, SpliceRegion), "splice_region should be SpliceRegion, not %s" % type(splice_region)

    plt.figure(
        figsize=[settings['width'], settings['height']]
    )

    plot_density(
        settings,                               # plot settings, untouched
        read_depths_dict=average_depths_dict,   # reads coverage
        splice_region=splice_region,            # Exon and transcript information
        show_gene=True,                         # decide whether display gene id in this plot
        # ordered_genotypes_list                # provide the allele information or label for subtitle in SplicePlot
    )

    logger.info("save to %s" % output_file_path)
    # plt.savefig(
    #     output_file_path,
    #     transparent=True,
    #     bbox_inches='tight'
    # )

    plt.show()