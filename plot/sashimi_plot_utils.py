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
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import pylab
from src.logger import logger
from src.SpliceRegion import SpliceRegion

from plot.density_plot_utils import plot_density_single
from plot.sideplot_utils import plot_sideplot
from plot.transcripts_plot_utils import plot_transcripts
from plot.utils import *


# Plot density for a series of bam files.
def plot_density(
        settings,
        read_depths_dict,
        splice_region,
        show_gene=False,
        title=None,
        no_bam=False,
        logtrans=None,
        distance_ratio=0.3,
        stack=False,
        show_side = True,
        side_strand_choice = None
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
    :param logtrans: y ticks log transformation or not, 2 -> log2; 10 -> log10，
                do not use `set_yscale` here, because there are junction under the x axis,
                and there coords do not convert into log by matplotlib, so it will cause a lot troubles
    :param distance_ratio: distance between transcript label and transcript line
    :param stack: whether to make stacked reads
    :param show_side: whether to show side plot
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
    reads = []

    if stack:   # modify information of reads
        for key, val in read_depths_dict.items():
            for i in val.reads:
                i.transcript = key.alias
                reads.append(i)

            exon_starts += val.exon_starts
            exon_ends += val.exon_ends

    # Get the right scalings
    graph_coords = get_scaling(
        tx_start, tx_end, strand,
        exon_starts, exon_ends,
        intron_scale, exon_scale,
        reverse_minus
    )

    # set the number of axis
    n_files = len(read_depths_dict) + ((len(transcripts) // 2) if len(transcripts) > 1 else 1)  + ((len(reads) // 4) if len(reads) > 1 else 1)

    number_of_sashimi = len(read_depths_dict)
    if show_side:
        n_files += len(read_depths_dict) * 2
        number_of_sashimi += len(read_depths_dict) * 2
    
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
    i = 0
    for sample_info in read_depths_dict.keys():
        average_read_depth = read_depths_dict[sample_info]

        curr_ax = plt.subplot(gs[i, :])

        if title is not None and i == 0:
            curr_ax.set_title(title, fontsize=10, loc="center")

        """
        Main body of sashimi
        @2019.01.07 
        """
        # Read sample label
        plot_density_single(
            read_depth_object=average_read_depth,
            graph_coords=graph_coords,
            ax_var=curr_ax,
            sample_info=sample_info,
            number_junctions=number_junctions,
            ny_ticks=nyticks,
            font_size=font_size,
            numbering_font_size=numbering_font_size,
            no_bam=no_bam,
            logtrans=logtrans,
            show_ylabel = show_ylabel,
            distance_between_label_axis = distance_between_label_axis,
            sites = splice_region.sites
        )

        """
        Draw side plot
        """
        if show_side:
            curr_ax = plt.subplot(gs[(i + 1):(i + 2), :])

            plot_sideplot(
                average_read_depth,
                sample_info,
                graph_coords,
                curr_ax,
                font_size=font_size,
                logtrans=logtrans,
                sites = splice_region.sites,
                strand_choice = side_strand_choice
            )

        i += 3 if show_side else 1
        if i == number_of_sashimi:
            # curr_ax = plt.subplot(gs[number_of_sashimi - 1, :])
            set_x_ticks(
                average_read_depth, curr_ax, 
                graph_coords, chromosome, strand, 
                logtrans = logtrans, nx_ticks = nxticks, 
                font_size=font_size
            )
        

   # Draw reads
    if len(reads) > 0:
        plt.subplot(gs[number_of_sashimi:len(reads) // 4 + len(read_depths_dict), :])
        plot_transcripts(
            tx_start=tx_start,
            transcripts=reads,
            graph_coords=graph_coords,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene,
            distance_ratio=distance_ratio,
            color="darkgrey"
        )
    # Draw gene structure
    """
    @2018.12.26
    add more subplots, based on the number of transcripts
    """
    if len(transcripts) > 0:
        plt.subplot(gs[number_of_sashimi + len(reads) // 4:, :]) # + 1 if splice_region.sequence else len(read_depths_dict)

        plot_transcripts(
            tx_start=tx_start,
            transcripts=transcripts,
            graph_coords=graph_coords,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene,
            distance_ratio=distance_ratio,
            ymax=len(transcripts) - .5,
            sites = splice_region.sites
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
        logtrans=None,
        distance_ratio=0.3,
        title=None,
        stack=False,
        show_side = False,
        side_strand_choice = None
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
    elif stack:
        num_reads = 0
        for val in average_depths_dict.values():
            num_reads += len(val.reads)

        temp_num = len(average_depths_dict) + \
            len(splice_region.transcripts) // 2 + \
            num_reads // 4

        height = settings['height'] * temp_num

        if splice_region.sequence:
            height += (settings["height"] * .2 if temp_num > 5 else .4)
    else:
        temp_num = (len(average_depths_dict) + len(splice_region.transcripts) // 2)

        height = settings['height'] * temp_num

        if splice_region.sequence:
            height += (settings["height"] * .2 if temp_num > 5 else .4)

    plt.figure(figsize=[settings['width'], height], dpi=dpi)

    plot_density(
        settings,                               # plot settings, untouched
        read_depths_dict=average_depths_dict,   # reads coverage
        splice_region=splice_region,            # Exon and transcript information
        show_gene=show_gene,                    # decide whether display gene id in this plot
        no_bam=no_bam,
        logtrans=logtrans,
        distance_ratio=distance_ratio,
        title=title,
        stack=stack,
        show_side = show_side,
        side_strand_choice = side_strand_choice
    )

    logger.info("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True,
        bbox_inches='tight'
    )

