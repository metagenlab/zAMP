"""
Command line interface for running zAMP
"""

import click
from snaketool_utils.cli_utils import OrderedCommands
from .utils import (
    get_version,
    snake_base,
    print_citation,
    common_options,
    db_options,
    run_options,
    run_snakemake,
)


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Snakemake pipeline designed for convenient, reproducible and scalable amplicon-based metagenomics

    For more options, run:
    zamp command --help"""
    pass


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@db_options
@common_options
def db(**kwargs):
    """
    Prepare database files for zAMP
    """
    merge_config = {"args": kwargs}

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base("DBprocess.Snakefile"),
        merge_config=merge_config,
        **kwargs,
    )


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@run_options
@common_options
def run(**kwargs):
    """
    Run zAMP
    """
    merge_config = {"args": kwargs}

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base("Snakefile"),
        merge_config=merge_config,
        **kwargs,
    )


@click.command()
def citation(**kwargs):
    """Print zAMP and tools citations"""
    print_citation()


cli.add_command(db)
cli.add_command(run)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
