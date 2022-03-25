#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang Yiming at 2021.03.16

This is migrated from Zhou's code
"""
import scipy.stats as sts

from plot.utils import *
from src.ReadDepth import ReadDepth
from src.SpliceRegion import SpliceRegion


def plot_sideplot(
        read_depth_object: ReadDepth,
        region: SpliceRegion,
        sample_info,
        ax_var,
        font_size=8,
        logtrans=None,
        strand_choice: str = None,
):
    """
    :return:
    """
    plus = read_depth_object.plus
    minus = read_depth_object.minus

    if logtrans == '2':
        plus = np.log2(plus + 1)
        minus = -np.log2(abs(minus) + 1)
    elif logtrans == '10':
        plus = np.log10(plus + 1)
        minus = -np.log10(abs(minus) + 1)
    else:
        pass

    maxheight = max(plus)
    minheight = min(minus)
    max_val = max(maxheight, abs(minheight))

    ymax = 1.1 * max_val
    ymin = 1.1 * -max_val

    for label, array_plot in zip(['plus', 'minus'], [plus, minus]):
        if strand_choice is not None and label != strand_choice:
            continue

        array_hist = np.repeat(region.graph_coords, np.abs(array_plot).astype(np.int))
        try:
            kde = sts.gaussian_kde(array_hist)
            fit_value = kde.pdf(region.graph_coords)
        except Exception:
            continue

        fit_value = fit_value / fit_value.max()
        if label == 'plus':
            ax_var.plot(region.graph_coords, fit_value * array_plot.max(), c=sample_info.color, lw=1)
            ax_var.bar(range(len(region)), array_plot, color=sample_info.color, rasterized=region.raster)
        else:
            ax_var.plot(region.graph_coords, fit_value * array_plot.min(), c=sample_info.color, lw=1)
            ax_var.bar(range(len(region)), array_plot, color=sample_info.color, rasterized=region.raster)

    # set the y limit
    # set y ticks, y label and label
    ax_var.set_ybound(lower=ymin, upper=ymax)
    universal_yticks = pylab.linspace(ymin, ymax, 3)

    curr_yticklabels = []
    for label in universal_yticks:
        curr_yticklabels.append("{}".format(int(label)))

    ax_var.set_yticks(universal_yticks)
    ax_var.set_yticklabels(curr_yticklabels, fontsize=font_size)

    ax_var.spines["left"].set_bounds(ymin, ymax)
    ax_var.yaxis.set_ticks_position('left')

    ax_var.spines["right"].set_visible(False)
    ax_var.spines['right'].set_color('none')
    ax_var.spines['top'].set_color('none')
    ax_var.spines['bottom'].set_color('none')

    pylab.xticks([])
    pylab.xlim(0, max(region.graph_coords))

    add_additional_background(region)

    print(region.raster)
    if region.raster:
        ax_var.set_rasterization_zorder(0)
    return ax_var


if __name__ == '__main__':
    pass
