#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.SpliceRegion import SpliceRegion
from src.logger import logger
from matplotlib import pylab


def plot_stroke(region: SpliceRegion, font_size: float = 5):
    strokes = sorted(region.stroke, key=lambda x: [x.start, x.end])

    for i, stroke in enumerate(strokes):
        pylab.hlines(y=i, xmin=stroke.start, xmax=stroke.end, color=stroke.color, lw=2)
        pylab.text(stroke.center - len(stroke.label) / 2, i + .2, stroke.label, fontsize=font_size)

    pylab.xlim(left=0, right=len(region))
    pylab.ylim(bottom=-1, top=len(strokes))
    pylab.box(on=False)
    pylab.xticks([])
    pylab.yticks([])


if __name__ == '__main__':
    pass
