#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
make sashimi plot
"""
import math

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

plt.switch_backend("TkAgg")


class Sashimi(object):
    u"""
    make the sashimi plot
    """

    def __init__(self, start, end, transcripts, coverages, sites, dpi=300):
        u"""
        init this class
        :param start:
        :param end:
        :param transcripts:
        """
        self.__x_start__ = 50
        self.__x_end__ = 150

        self.start = int(start)
        self.end = int(end)
        self.scale_size = (self.end - self.start) / (self.__x_end__ - self.__x_start__)
        self.sites = [int(x) for x in sites.split(",")] if sites else None

        self.transcripts = [x for x in transcripts] if transcripts is not None else []

        self.coverages = coverages

        self.colors = sns.color_palette("husl", 8)

        self.width = self.__x_end__ + 10
        self.height = len(self.transcripts) * 10 + len(self.coverages) * 60 + 60

        # set fig size and axises
        self.fig, self.ax = plt.subplots(dpi=dpi)
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(0, self.height)
        self.ax.set_axis_off()

        # plot boundary
        self.boundary = len(self.transcripts) * 10 + 60

        # plot figure
        self.__plot_border__()
        self.__plot_coverage__()
        self.__plot_transcripts__()

    def __plot_border__(self):
        u"""
        plot borders of reads coverage
        :return:
        """
        self.ax.plot(
            [self.__x_start__, self.__x_end__],
            [self.boundary - 30, self.boundary - 30],
            linestyle="-",
            color="black",
            linewidth=0.5
        )

        for i in range(len(self.coverages)):
            bottom = self.boundary + 60 * i
            # self.ax.plot(
            #     [self.__x_start__, self.__x_end__],
            #     [bottom, bottom],
            #     linestyle="-",
            #     color="black",
            #     linewidth=0.5
            # )

            self.ax.plot(
                [self.__x_start__, self.__x_start__],
                [bottom, bottom + 50],
                linestyle="-",
                color="black",
                linewidth=0.5
            )

            # self.ax.plot(
            #     [self.__x_end__, self.__x_end__],
            #     [bottom, bottom + 50],
            #     linestyle="-",
            #     color="black",
            #     linewidth=0.5
            # )

    def __scale_exons_sites__(self, exons):
        u"""
        scale exons site to 0-100
        :param exons:
        :return:
        """
        if exons[0] < self.start:
            exons[0] = self.start

        res = []
        for i in exons:
            tmp = (i - self.start) / self.scale_size + self.__x_start__

            if tmp > self.__x_end__:
                break

            res.append(tmp)

        return res

    def __plot_transcripts__(self):
        u"""
        plot the transcripts
        :return:
        """
        for i in range(len(self.transcripts)):
            self.ax.text(
                5,
                i * 10 + 4,
                self.transcripts[i][0],
                fontsize=5
            )

            # plot each transcript
            self.ax.plot(
                [self.__x_start__, self.__x_end__],
                [i * 10 + 5, i * 10 + 5],
                linestyle="-",
                linewidth=0.5,
                color="black",
                marker=9 if self.transcripts[i][1] == "+" else 8,
                markersize=2,
                markevery=0.1
            )
            # add exons
            exons = self.__scale_exons_sites__(self.transcripts[i][2])

            for j in range(0, len(exons) - 1, 2):
                rect = plt.Rectangle(
                    (exons[j], i * 10 + 3),      # bottom and left
                    exons[j + 1] - exons[j],
                    5,
                    color="black"
                )
                self.ax.add_patch(rect)

    def __plot_junction__(self, junction, count, y, color, bottom, up=True):
        u"""
        plot junctions
        :param junction: [start, end]
        :param count: number of this junctions
        :param y: y axises, calculated by reads coverages
        :param color: color of this junction line
        :param up: Boolean -> plot junction on top of coverage
        :return:
        """

        # get the x and y axis
        start = (junction[0] - self.start) / self.scale_size + self.__x_start__
        end = (junction[1] - self.start) / self.scale_size + 1 + self.__x_start__
        middle = (start + end) / 2

        # get the x and y axis
        x_axis = np.linspace(0, np.pi, 201)

        if up:
            bottom = min([y[junction[0] - self.start - 1], y[junction[1] - self.start + 1]])
            y_axis = np.sin(x_axis) * 5 + bottom
        else:
            y_axis = np.sin(x_axis) * -5 + bottom

        # print(y_axis)
        self.ax.plot(
            x_axis * (end - start) / np.pi + start,
            y_axis,
            c=color,
            linewidth=0.2
        )

        # get text position
        self.ax.text(
            middle,
            max(y_axis) + 1 if up else min(y_axis) - 2,
            str(count),
            fontsize=3
        )

    def __plot_coverage__(self):
        u"""
        plot the reads coverages
        """
        for i in range(len(self.coverages)):
            coverage = self.coverages[i].coverage
            junctions = self.coverages[i].introns

            scale_y = max(coverage) / 40
            bottom = self.boundary + 60 * i + 5

            x, y1, y2 = [], [], []
            for j in range(len(coverage)):

                x.append(j / self.scale_size + self.__x_start__)
                if scale_y == 0:
                    y1.append(bottom)
                else:
                    y1.append(coverage[j] / scale_y + bottom)
                y2.append(bottom)
                # y = self.boundary + 60 * i

            self.ax.fill_between(
                x,
                y1,
                y2,
                facecolor=self.colors[i]
            )

            # plot annotation lines
            self.__plot_sites__(y2[0])

            # plot junctions
            up = True
            for junction, count in junctions.items():
                self.__plot_junction__(
                    junction=junction,
                    count=count,
                    y=y1,
                    color=self.colors[i],
                    bottom=bottom,
                    up=up
                )
                up = not up

            self.ax.text(
                -15,
                bottom + 25,
                self.coverages[i].name,
                fontsize=8
            )

        pass

    def __plot_sites__(self, y):
        u"""
        add annotation line
        """
        if not self.sites:
            return
        for i in self.sites:
            i = (i - self.start) / self.scale_size + self.__x_start__
            self.ax.plot(
                [i, i],
                [y, y + 60],
                c="grey",
                linestyle="--"
            )

    def save(self, outfile):
        u"""
        :param outfile:
        :return:
        """
        if outfile is None:
            plt.show()
        else:
            self.fig.savefig(outfile)
            plt.close(self.fig)

