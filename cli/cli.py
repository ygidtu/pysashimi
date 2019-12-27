VERSION = "1.4.0"
LABEL = "pySashimi"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

import click

from cli.no_bam import no_bam
from cli.normal import normal
from cli.pipeline import pipeline
from cli.line import line

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


cli.add_command(normal)
cli.add_command(no_bam)
cli.add_command(pipeline)
cli.add_command(line)