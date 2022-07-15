#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""


class AxLabel(object):
    def __init__(self, ax, label):
        self.Ax = ax
        self.Label = label

    def __hash__(self):
        return hash(self.Ax)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


if __name__ == '__main__':
    pass
