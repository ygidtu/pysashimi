#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import seaborn as sns

from src.Bigwig import Bigwig


def plot_heatmap(bigwig: Bigwig, ax_var, cbar_ax, font_size: float = 4, distance_between_label_axis: float = 0):
    u"""
    :return:
    """
    if bigwig.data is None:
        raise ValueError("Please prepare Bigwig first")

    sns.heatmap(
        bigwig.data,
        ax=ax_var, cmap=bigwig.color,
        cbar_ax=cbar_ax,
        xticklabels=False,
        yticklabels=False,
        center=True,
        rasterized=bigwig.raster
    )

    ax_var.tick_params(axis='both', which='major', labelsize=font_size)
    cbar_ax.tick_params(labelsize=font_size)

    if bigwig.alias:
        ax_var.set_ylabel(
            bigwig.alias,
            fontsize=font_size,
            va="center",
            labelpad=distance_between_label_axis,  # the distance between y label with axis
            rotation="horizontal"
        )

    if bigwig.raster:
        ax_var.set_rasterization_zorder(0)


if __name__ == '__main__':
    pass
