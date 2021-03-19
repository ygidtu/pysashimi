#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang at 2021.03.16

Seperate the  single density plot from sashimi_plot_utils
"""

import numpy
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from plot.utils import *


def plot_density_single(
        read_depth_object,
        chromosome,
        strand,
        graph_coords,
        ax_var,
        sample_info,
        distance_between_label_axis,
        number_junctions=True,
        show_ylabel = True,
        nx_ticks=4,
        ny_ticks=4,
        font_size=6,
        numbering_font_size=6,
        no_bam=False,
        logtrans=None,
        sites = None
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
    :param sample_info:
    :param number_junctions:
    :param nx_ticks:
    :param font_size:
    :param numbering_font_size:
    :param no_bam:
    :param logtrans:
    :return:
    """
    # Round up
    max_used_y_val = read_depth_object.max
    fake_y_min = - 0.5 * max_used_y_val
    """
    @2018.12.20 if max_used_y_val is odd, plus one, for better look
    @2019.01.07 flush universal_y_ticks if not share_y
    """
    if max_used_y_val % 2 == 1:
        max_used_y_val += 1

    # @2018.12.19
    # tx_end = read_depth_object.high
    # chrom = read_depth_object.chrm

    wiggle = read_depth_object.wiggle

    jxns = read_depth_object.junctions_dict

    y_max = max_used_y_val
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
            color=sample_info.color,
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
        ss1_idx, ss1_modified = get_limited_index(leftss - read_depth_object.start, len(graph_coords))
        ss2_idx, ss2_modified = get_limited_index(rightss - read_depth_object.start, len(graph_coords))

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
            a, ec=sample_info.color,
            lw=line_width + 0.2, fc='none'
        )

        ax_var.add_patch(p)

    # set y ticks, y label and label

    ax_var.set_ybound(lower=fake_y_min, upper=1.2 * max_used_y_val)
    ax_var.spines["left"].set_bounds(0, max_used_y_val)
    ax_var.spines["right"].set_color('none')

    universal_y_ticks = pylab.linspace(0, max_used_y_val, ny_ticks + 1)

    set_y_ticks(
        ax_var, sample_info, 
        universal_y_ticks, 
        distance_between_label_axis,
        font_size = font_size, 
        logtrans = logtrans, 
        show_ylabel = show_ylabel, 
        no_bam=no_bam
    )
    set_label(ax_var, sample_info, graph_coords, max_used_y_val, font_size=font_size)
    set_indicator_lines(read_depth_object, ax_var, graph_coords, sites, max_used_y_val)


    # Format plot
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')

    ax_var.spines['bottom'].set_color('none')
    pylab.xticks([])

    pylab.xlim(0, max(graph_coords))
