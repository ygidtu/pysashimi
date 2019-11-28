#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2018.12.16

Main function to plot sashimi plot
"""
import click

from cli.no_bam import no_bam
from cli.normal import normal
from cli.pipeline import pipeline
from cli.line import line


VERSION = "1.4.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
)
@click.version_option(VERSION, message="Current version %(version)s")
def main():
    u"""
    Welcome

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :return:
    """

    pass


if __name__ == '__main__':
    main.add_command(normal)
    main.add_command(no_bam)
    main.add_command(pipeline)
    main.add_command(line)
    main()
    pass


