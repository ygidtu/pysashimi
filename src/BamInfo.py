#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
from typing import Dict, List, Optional


class BamInfo(object):
    def __init__(self, alias, title, label, path, color, barcodes = None):
        self.alias = alias
        self.title = title
        self.label = label
        self.path = [path]
        self.color = color
        self.barcodes = self.set_barcodes(barcodes)

    def set_barcodes(self, barcodes: Optional[List[str]]) -> Dict:
        u"""
        seperate barcodes by it's first character to reduce set size
        :params barcodes: list or set of barcodes
        """
        res = {}

        if barcodes is not None:
            for b in barcodes:
                if b:
                    f = b[:min(3, len(b))]

                    if f not in res.keys():
                        res[f] = set()

                    res[f].add(b)

        return res

    def has_barcode(self, barcode: str) -> bool:
        u"""
        check whether contains barcodes
        :param barcode: barcode string
        """
        if barcode:
            f = barcode[:min(3, len(barcode))]

            temp = self.barcodes.get(f, set())

            return barcode in temp
        return False

    def empty_barcode(self) -> bool:
        u"""
        check whether this bam do not contains any barcodes
        
        """
        count = 0

        for i in self.barcodes.values():
            count += len(i)

            if count > 0:
                return False

        return True


    def __hash__(self):
        return hash(self.alias)

    def __str__(self) -> str:

        temp = []

        for x in [self.alias, self.title, self.label, self.path, self.color]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return "\t".join(temp)

    def __eq__(self, other) -> bool:
        return self.__hash__() == other.__hash__()

    def to_csv(self) -> str:
        temp = []

        for x in [self.alias, self.title, self.label, self.path, self.color]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return ",".join(temp)

    def __add__(self, other):
        self.path += other.path
        
        for i, j in other.barcodes.items():
            if i not in self.barcodes.keys():
                self.barcodes[i] = j
            else:
                self.barcodes[i] |= j

        return self


if __name__ == '__main__':
    pass
