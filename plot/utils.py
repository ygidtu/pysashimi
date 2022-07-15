#!/usr/bin/env python3
# -*- utf-8 -*-
u"""
Created at 2021.03.16

basic func
"""
import numpy as np
from matplotlib import pylab

from src.ReadDepth import ReadDepth
from src.SpliceRegion import SpliceRegion
from src.logger import logger


def get_limited_index(num, length):
    u"""
    Created by ygidtu at 2018.12.19

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
    u"""
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = pylab.array(p0)
    p1 = pylab.array(p1)
    p2 = pylab.array(p2)
    p3 = pylab.array(p3)
    return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + 3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3


def set_x_ticks(read_depth_object: ReadDepth, ax_var, region: SpliceRegion, logtrans=None,
                nx_ticks=4, font_size=4):
    ax_var.xaxis.set_ticks_position('bottom')

    # @2018.12.19 unnecessary text in figure

    xlabel = f'Genomic coordinate ({region.chromosome}), "{region.strand}" strand'

    if logtrans in (2, 10):
        xlabel = f"{xlabel}, y axis is log{logtrans} transformed"

    pylab.xlabel(xlabel, fontsize=font_size)

    bk = 1
    if not read_depth_object.sequence:
        bk = len(region) // nx_ticks

    linspace, ticks = [], []
    for i in range(0, len(region), bk):
        linspace.append(region.graph_coords[i])

        temp_txs = read_depth_object.start + i
        if read_depth_object.sequence:
            if (i - read_depth_object.start) % nx_ticks == 0:
                temp_txs = "{}\n{}".format(read_depth_object.start + i, read_depth_object.sequence[i])
            else:
                temp_txs = "\n{}".format(read_depth_object.sequence[i])
        ticks.append(temp_txs)

    # ax_var.spines['right'].set_color('none')
    # ax_var.spines['top'].set_color('none')
    # pylab.yticks([])
    ax_var.spines['bottom'].set_color('black')
    ax_var.xaxis.set_ticks_position('bottom')
    pylab.xticks(linspace, ticks, fontsize=font_size)


def set_y_ticks(ax_var, sample_info, universal_y_ticks, distance_between_label_axis, font_size=5, logtrans=None,
                show_ylabel=True, no_bam=False):
    u"""
    The y ticks are formatted here
    @2019.03.31 add little check here to make sure the y axis shows the real value
    """
    curr_y_tick_labels = []

    for label in universal_y_ticks:
        if logtrans in (2, 10):
            label = np.power(logtrans, label)

        if label <= 0:
            # Exclude label for 0
            curr_y_tick_labels.append("")
        else:
            curr_y_tick_labels.append("%.1f" % label if label % 1 != 0 else "%d" % label)

    if not no_bam:
        ax_var.set_yticks(universal_y_ticks)
        ax_var.set_yticklabels(
            curr_y_tick_labels,
            fontsize=font_size
        )

        ax_var.yaxis.set_ticks_position('left')
    else:
        u"""
        @2019.01.04

        If there is no bam file, draw a blank y-axis 
        """
        ax_var.set_yticks([])
        ax_var.yaxis.set_ticks_position("none")

    """
    Plot y labels
    @2018.12.20 using BAM label as ylabel

    @2019.01.04 change the standards of distance between ylabel and y-axis
    """
    if show_ylabel:
        ax_var.set_ylabel(
            sample_info.alias,
            fontsize=font_size,
            va="center",
            labelpad=distance_between_label_axis,  # the distance between ylabel with axis
            rotation="horizontal"
        )


def set_label(ax_var, sample_info, graph_coords, ymax, font_size=5):
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

    t = ax_var.text(
        max(graph_coords),
        ymax,
        curr_label,
        fontsize=font_size,
        va='bottom',
        ha='right',
        color=sample_info.color
    )

    # @218.12.19 set transparent background
    t.set_bbox(dict(alpha=0))


def add_additional_background(region: SpliceRegion):
    if region.sites is None:
        return
    _, _, y1, y2 = pylab.axis()

    for site, color in region.sites.items():
        try:
            pylab.vlines(
                x=site,
                ymin=y1,
                ymax=y2,
                color=color,
                linestyles="dashed",
                lw=0.5
            )
        except IndexError as err:
            logger.warning("Indicator line is out of bound: " + str(err))

    if region.focus is not None:
        for left, right in region.focus.items():
            try:
                fill_x = [left, right, right, left]

                _, _, y1, y2 = pylab.axis()
                fill_y = [y1, y1, y2, y2]

                pylab.fill(fill_x, fill_y, alpha=0.1, color='grey')
            except IndexError as err:
                logger.warning("focus region is out of bound: " + str(err))


if __name__ == '__main__':
    pass
