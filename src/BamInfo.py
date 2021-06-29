#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""


class BamInfo(object):
    def __init__(self, alias, title, label, path, color, barcodes = None):
        self.alias = alias
        self.title = title
        self.label = label
        self.path = [path]
        self.color = color
        self.barcodes = barcodes

        if self.barcodes is None:
            self.barcodes = []

    def __hash__(self):
        return hash(self.alias)

    def __str__(self):

        temp = []

        for x in [self.alias, self.title, self.label, self.path, self.color]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return "\t".join(temp)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def to_csv(self):
        temp = []

        for x in [self.alias, self.title, self.label, self.path, self.color]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return ",".join(temp)

    def __add__(self, other):
        self.path += other.path
        temp = list(self.barcodes)
        temp += list(other.barcodes)
        self.barcodes = temp

        return self


if __name__ == '__main__':
    pass
