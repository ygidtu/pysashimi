#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.11.20

Make line
"""
from typing import Dict

from matplotlib import pylab

from src.BamInfo import BamInfo
from src.ReadDepth import ReadDepth
from src.SpliceRegion import SpliceRegion


def plot_line_single(
        read_depth_dict: Dict[BamInfo, ReadDepth],
        region: SpliceRegion,
        ax_var,
        show_x_axis=True,
        nx_ticks=4,
        font_size=6,
        logtrans=None
):
    u"""
    @2018.12.19 remove unnecessary x label
    @2018.12.19 replace junc_comp_function with a new Junction class,
                due to cmp function is removed in Python3
    """

    # Plot all reads depth as line
    for bam, obj in read_depth_dict.items():
        wiggle = obj.wiggle

        # Reduce memory footprint by using incremented graph coords.
        compressed_x = []
        compressed_wiggle = []

        for i in range(len(region)):
            compressed_wiggle.append(wiggle[i])
            compressed_x.append(region.graph_coords[i])

        pylab.plot(compressed_x, compressed_wiggle, color=bam.color, lw=1)

    # Format plot
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')

    if show_x_axis:
        ax_var.xaxis.set_ticks_position('bottom')

        # @2018.12.19 unnecessary text in figure

        pylab.xlabel(
            region.x_label(logtrans),
            fontsize=font_size
        )

        bk = len(region) // nx_ticks

        linspace, ticks = [], []
        for i in range(0, len(region), bk):
            linspace.append(region.graph_coords[i])
            ticks.append(region.start + i)

        pylab.xticks(linspace, ticks, fontsize=font_size)

    else:
        ax_var.spines['bottom'].set_color('none')
        pylab.xticks([])

    pylab.xlim(0, max(region.graph_coords))


if __name__ == '__main__':
    pass
