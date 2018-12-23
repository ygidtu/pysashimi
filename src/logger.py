#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by ygidtu at 2018.12.21

logger section
"""
import logging


def set_logging(log_name):
    u"""
    set logging handler
    Created by Zhang yiming at 2018.11.13
    :return:
    """
    sh = logging.StreamHandler()

    # [%(module)s:%(lineno)d]
    formatter = logging.Formatter(
        fmt="[%(asctime)s] - [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)

    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    log.addHandler(sh)
    return log


logger = set_logging("pySplice")

if __name__ == '__main__':
    pass
