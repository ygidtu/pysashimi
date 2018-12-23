#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import matplotlib.pyplot as plt

if __name__ == '__main__':

    t = plt.text(0, 0, "test", alpha=0.5)
    t.set_bbox(dict(alpha=0))
    plt.show()
