#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
make sashimi plot
"""
import math

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

plt.switch_backend("Agg")


class Sashimi(object):
    u"""
    make the sashimi plot
    """

    def __init__(self, start, end, transcripts, coverages):
        u"""
        init this class
        :param start:
        :param end:
        :param transcripts:
        """

        self.start = int(start)
        self.end = int(end)
        self.scale_size = (self.end - self.start) / 100

        self.transcripts = [x for x in transcripts] if transcripts is not None else []

        self.coverages = coverages

        self.colors = sns.color_palette("husl", 8)

        # left bottom anchor of sashimi plot
        self.__x_start__ = 50
        self.__x_end__ = 150

        self.width = self.__x_end__ + 10
        self.height = len(self.transcripts) * 10 + len(self.coverages) * 60

        # set fig size and axises
        self.fig, self.ax = plt.subplots(dpi=300)
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(0, self.height)
        self.ax.set_axis_off()

        # plot boundary
        self.boundary = len(self.transcripts) * 10 + 10

        # plot figure
        self.__plot_border__()
        self.__plot_coverage__()
        self.__plot_transcripts__()

    def __plot_border__(self):
        u"""
        plot borders of reads coverage
        :return:
        """
        for i in range(len(self.coverages)):
            bottom = self.boundary + 60 * i
            self.ax.plot(
                [self.__x_start__, self.__x_end__],
                [bottom, bottom],
                linestyle="-",
                color="black",
                linewidth=0.5
            )

            self.ax.plot(
                [self.__x_start__, self.__x_start__],
                [bottom, bottom + 50],
                linestyle="-",
                color="black",
                linewidth=0.5
            )

            self.ax.plot(
                [self.__x_end__, self.__x_end__],
                [bottom, bottom + 50],
                linestyle="-",
                color="black",
                linewidth=0.5
            )

    def __scale_exons_sites__(self, exons):
        u"""
        scale exons site to 0-100
        :param exons:
        :return:
        """
        if exons[0] < self.start:
            exons[0] = self.start
        if exons[-1] > self.end:
            exons[-1] = self.end

        return [(x - self.start) / self.scale_size + self.__x_start__ for x in exons]

    def __plot_transcripts__(self):
        u"""
        plot the transcripts
        :return:
        """
        for i in range(len(self.transcripts)):
            self.ax.text(
                -15,
                i * 10 + 4,
                self.transcripts[i][0],
                fontsize=8
            )

            # plot each transcript
            self.ax.plot(
                [self.__x_start__, self.__x_end__],
                [i * 10 + 5, i * 10 + 5],
                linestyle="-",
                linewidth=2,
                color="black",
                marker=9 if self.transcripts[i][1] == "+" else 8,
                markevery=0.1
            )
            # add exons
            exons = self.__scale_exons_sites__(self.transcripts[i][2])

            for j in range(0, len(exons) - 1, 2):
                rect = plt.Rectangle(
                    [exons[j], i * 10 + 3],
                    exons[j + 1] - exons[j],
                    5,
                    color="black"
                )
                self.ax.add_patch(rect)

    def __plot_junction__(self, junction, count, y, color):
        u"""
        plot junctions
        :param junction: [start, end]
        :param count: number of this junctions
        :param y: y axises, calculated by reads coverages
        :param color: color of this junction line
        :return:
        """
        bottom = min([y[junction[0] - self.start - 1], y[junction[1] - self.start + 1]])

        # get the x and y axis
        start = (junction[0] - self.start) / self.scale_size + self.__x_start__
        end = (junction[1] - self.start) / self.scale_size + 1 + self.__x_start__
        middle = (start + end) / 2

        # get the x and y axis
        x_axis = np.linspace(0, np.pi, 201)
        y_axis = np.sin(x_axis) * 5 + bottom

        # print(y_axis)
        self.ax.plot(
            x_axis * (end - start) / np.pi + start,
            y_axis,
            c=color,
            linewidth=math.log10(count + 1)
        )

        # get text position
        self.ax.text(middle, max(y_axis) + 1, str(count), fontsize=5)

    def __plot_coverage__(self):
        u"""
        plot the reads coverages
        """
        for i in range(len(self.coverages)):
            coverage = self.coverages[i].coverage
            junctions = self.coverages[i].introns

            scale_y = max(coverage) / 40
            bottom = self.boundary + 60 * i

            x, y1, y2 = [], [], []
            for j in range(len(coverage)):

                x.append(j / self.scale_size + self.__x_start__)
                y1.append(coverage[j] / scale_y + bottom)
                y2.append(self.boundary + 60 * i)
                # y = self.boundary + 60 * i

            self.ax.fill_between(
                x,
                y1,
                y2,
                facecolor=self.colors[i]
            )

            for junction, count in junctions.items():
                self.__plot_junction__(
                    junction=junction,
                    count=count,
                    y=y1,
                    base=bottom,
                    color=self.colors[i]
                )

            self.ax.text(
                -15,
                bottom + 25,
                self.coverages[i].name,
                fontsize=8
            )

        pass

    def save(self, outfile):
        u"""
        :param outfile:
        :return:
        """
        self.fig.savefig(outfile)
        plt.close(self.fig)

