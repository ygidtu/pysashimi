#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.12.27
"""
import click
from src.logger import init_logger

from cli.line import line
from cli.plot import plot

VERSION = "1.5.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

init_logger("INFO")

@click.group(
    context_settings=CONTEXT_SETTINGS,
)
@click.version_option(VERSION, message="Current version %(version)s")
def cli():
    u"""
    Welcome

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :return:
    """

    pass


cli.add_command(plot)
cli.add_command(line)
