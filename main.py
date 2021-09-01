#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

cli function to plot sashimi plot
"""
import matplotlib as mpl
import matplotlib.font_manager

mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42

if any(["Arial" in f.name for f in matplotlib.font_manager.fontManager.ttflist ]):
    mpl.rcParams['font.family'] = 'Arial'

from cli.cli import cli


if __name__ == '__main__':
    cli()
    pass


