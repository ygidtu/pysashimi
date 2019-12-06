#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""


class BamInfo(object):
    def __init__(self, alias, title, label, path, color):
        self.alias = alias
        self.title = title
        self.label = label
        self.path = path
        self.color = color

    def __hash__(self):
        return hash(self.path)

    def __str__(self):

        temp = []

        for x in [self.alias, self.title, self.label, self.path, self.color]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return "\t".join(temp)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()