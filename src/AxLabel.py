#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""


class AxLabel(object):
    def __init__(self, Ax, Label):
        self.Ax = Ax
        self.Label = Label

    def __hash__(self):
        return hash(self.Ax)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
