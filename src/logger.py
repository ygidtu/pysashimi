#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.11 by Zhang
"""

import logging

from rich.logging import RichHandler

logger = logging.getLogger("rich")


def init_logger(level="NOTSET"):
    logging.basicConfig(
        level=level, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
    )
    global logger
    logger = logging.getLogger("rich")


if __name__ == '__main__':
    pass
