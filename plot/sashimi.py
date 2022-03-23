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
9. add parameter to decide if plot by shared y axis
10. fix transcripts display issues
"""
from typing import Dict

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from plot.density import plot_density_single
from plot.heatmap import plot_heatmap
from plot.line import plot_line_single
from plot.sideplot import plot_sideplot
from plot.stroke import plot_stroke
from plot.transcripts import plot_transcripts
from plot.utils import *
from src.BamInfo import BamInfo
from src.ReadDepth import ReadDepth
from src.SpliceRegion import SpliceRegion


# Plot density for a series of bam files.
def plot(
        settings,
        read_depths_dict: Dict[BamInfo, ReadDepth],
        splice_region: SpliceRegion,
        show_gene=False,
        title=None,
        no_bam=False,
        logtrans=None,
        distance_ratio=0.3,
        stack=False,
        show_side=True,
        side_strand_choice=None,
        draw_line: bool = False,
        **kwargs
):
    u"""
    Several modifications were taken

    1. Due to the change of mRNAObjects, the plot_mRNAs need to be modified
    2. ordered_genotypes_list used to plot the different allele specific situation,
        for now, I changed this to plot multiple BAM files

    :param draw_line: draw multiple lines in single plot
    :param side_strand_choice:
    :param settings:
    :param read_depths_dict:
    :param splice_region:
    :param show_gene: Boolean value to decide whether show gene id in this graph,
                    used by plot_transcripts()
    :param title: the title of this plot
    :param no_bam:
    :param logtrans: y ticks log transformation or not, 2 -> log2; 10 -> log10，
                do not use `set_yscale` here, because there are junction under the x axis,
                and there coords do not convert into log by matplotlib, so it will cause lot troubles
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
    transcripts = splice_region.transcripts
    reads = []

    splice_region.init_graph_coords(read_depths_dict if stack else None, intron_scale, exon_scale, reverse_minus)

    if stack:  # modify information of reads
        for key, val in read_depths_dict.items():
            for i in val.reads:
                i.transcript = key.alias
                reads.append(i)

    # set the number of axis
    number_of_sashimi = 1 if draw_line else len(read_depths_dict)
    if kwargs.get("bigwigs"):
        number_of_sashimi += len(kwargs.get("bigwigs"))

    n_files = number_of_sashimi

    if len(transcripts) > 1:
        n_files += (len(transcripts) // 2)

    if stack:
        n_files += (len(reads) // 4) if len(reads) > 1 else 1

    if show_side and not draw_line:
        n_files += len(read_depths_dict) * 2
        number_of_sashimi += len(read_depths_dict) * 2

    if splice_region.stroke:
        n_files += 1

    gs = gridspec.GridSpec(n_files, 1, hspace=.15, wspace=.7)
    if kwargs.get("bigwigs"):
        gs = gridspec.GridSpec(n_files, 2, width_ratios=(.99, .01), wspace=0.01, hspace=.15)

    """
    @2019.01.07
    calculate the distance between ylabel and y axis
    """
    distance_between_label_axis = max([len(x) if isinstance(
        x, str) else len(x.alias) for x in read_depths_dict.keys()]) * 2.5

    for i, j in enumerate(kwargs.get("bigwigs", [])):
        plot_heatmap(
            j,
            ax_var=plt.subplot(gs[i, 0]),
            cbar_ax=plt.subplot(gs[i, 1]),
            font_size=font_size,
            distance_between_label_axis=distance_between_label_axis
        )

    u"""
    @2018.12.19

    This part of code, used to plot different allele specific, but I this to plot multiple BAM files
    """
    i = len(kwargs.get("bigwigs", []))

    for sample_info in read_depths_dict.keys():
        average_read_depth = read_depths_dict[sample_info]

        curr_ax = plt.subplot(gs[i, 0])

        if title is not None and i == 0:
            curr_ax.set_title(title, fontsize=10, loc="center")

        """
        Main body of sashimi
        @2019.01.07 
        """
        # Read sample label
        if draw_line:
            plot_line_single(
                read_depth_dict=read_depths_dict,
                ax_var=curr_ax,
                font_size=font_size,
                logtrans=logtrans,
                region=splice_region
            )
        else:
            plot_density_single(
                read_depth_object=average_read_depth,
                ax_var=curr_ax,
                sample_info=sample_info,
                number_junctions=number_junctions,
                ny_ticks=nyticks,
                font_size=font_size,
                numbering_font_size=numbering_font_size,
                no_bam=no_bam,
                logtrans=logtrans,
                show_ylabel=show_ylabel,
                distance_between_label_axis=distance_between_label_axis,
                region=splice_region
            )

        """
        Draw side plot
        """
        if show_side and not draw_line:
            curr_ax = plt.subplot(gs[(i + 1):(i + 2), 0])

            plot_sideplot(
                average_read_depth,
                region=splice_region,
                sample_info=sample_info,
                ax_var=curr_ax,
                font_size=font_size,
                logtrans=logtrans,
                strand_choice=side_strand_choice
            )

        i += 3 if show_side else 1
        if i == number_of_sashimi:
            set_x_ticks(
                average_read_depth, curr_ax,
                region=splice_region,
                logtrans=logtrans, nx_ticks=nxticks,
                font_size=font_size
            )

    if splice_region.sequence:
        i += 1

    # Draw reads
    if len(reads) > 0:
        plt.subplot(gs[i:len(reads) // 4, 0])
        plot_transcripts(
            transcripts=reads,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene,
            distance_ratio=distance_ratio,
            color="darkgrey",
            region=splice_region
        )
        i += len(reads) // 4

    # Draw gene structure
    """
    @2018.12.26
    add more subplots, based on the number of transcripts
    """
    if len(transcripts) > 0:
        # + 1 if splice_region.sequence else len(read_depths_dict)
        plt.subplot(gs[i:(i + len(transcripts) // 2), 0])

        plot_transcripts(
            transcripts=transcripts,
            reverse_minus=reverse_minus,
            font_size=font_size,
            show_gene=show_gene,
            distance_ratio=distance_ratio,
            region=splice_region
        )
        i += len(transcripts) // 2

    if splice_region.stroke:
        plt.subplot(gs[i:, 0])
        plot_stroke(splice_region, font_size=font_size)


def save_fig(
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
        show_side=False,
        side_strand_choice=None,
        draw_line: bool = False,
        **kwargs
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

    assert isinstance(
        splice_region, SpliceRegion), "splice_region should be SpliceRegion, not %s" % type(splice_region)

    u"""
    @2019.01.04
    If there is no bam, reduce the height of figure
    """
    wigs_num = len(kwargs.get("bigwigs", [])) * 2

    if no_bam:
        height = settings['height'] * \
                 (len(average_depths_dict) + len(splice_region.transcripts) + wigs_num) // 2
    elif stack:
        num_reads = 0
        for val in average_depths_dict.values():
            num_reads += len(val.reads)

        temp_num = len(average_depths_dict) + len(splice_region.transcripts) // 2 + wigs_num + num_reads // 4

        height = settings['height'] * temp_num

        if splice_region.sequence:
            height += (settings["height"] * .2 if temp_num > 5 else .4)
    else:
        temp_num = (len(average_depths_dict) +
                    len(splice_region.transcripts) // 2)

        height = settings['height'] * temp_num + wigs_num

        if splice_region.sequence:
            height += (settings["height"] * .2 if temp_num > 5 else .4)

    plt.figure(figsize=[settings['width'], height], dpi=dpi)

    plot(
        settings,  # plot settings, untouched
        read_depths_dict=average_depths_dict,  # reads coverage
        splice_region=splice_region,  # Exon and transcript information
        show_gene=show_gene,  # decide whether display gene id in this plot
        no_bam=no_bam,
        logtrans=logtrans,
        distance_ratio=distance_ratio,
        title=title,
        stack=stack,
        show_side=show_side,
        side_strand_choice=side_strand_choice,
        draw_line=draw_line,
        **kwargs
    )

    logger.info("save to %s" % output_file_path)
    plt.savefig(
        output_file_path,
        transparent=True,
        bbox_inches='tight'
    )


if __name__ == '__main__':
    pass
