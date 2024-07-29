import os
import click
import yaml
import sys
import subprocess
import yaml
import json

from snaketool_utils.cli_utils import (
    msg,
    read_config,
    update_config,
    copy_config,
    echo_click,
    msg_box,
)


def validate_dict(ctx, param, value):
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        raise click.BadParameter(
            "Input should be a valid JSON string representing a dictionary."
        )


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("zamp.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("zamp.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def run_snakemake(
    configfile=None,
    system_config=None,
    snakefile_path=None,
    merge_config=None,
    threads=1,
    use_singularity=False,
    singularity_prefix=None,
    use_conda=False,
    conda_prefix=None,
    snake_default=None,
    snake_args=[],
    profile=None,
    workflow_profile=None,
    system_workflow_profile=None,
    log=None,
    **kwargs,
):
    """Run a Snakefile!

    Args:
        configfile (str): Filepath of config file to pass with --configfile
        system_config (str): Filepath of system config to copy if configfile not present
        snakefile_path (str): Filepath of Snakefile
        merge_config (dict): Config values to merge with your config file
        threads (int): Number of local threads to request
        use_conda (bool): Snakemake's --use-conda
        conda_prefix (str): Filepath for Snakemake's --conda-prefix
        snake_default (list): Snakemake args to pass to Snakemake
        snake_args (list): Additional args to pass to Snakemake
        profile (str): Name of Snakemake profile
        workflow_profile (str): Name of Snakemake workflow-profile
        system_workflow_profile (str): Filepath of system workflow-profile config.yaml to copy if not present
        log (str): Log file for writing STDERR
        **kwargs:

    Returns (int): Exit code
    """

    snake_command = ["snakemake", "-s", snakefile_path]

    # if using a configfile
    if configfile:
        # copy sys default config if not present
        copy_config(configfile, system_config=system_config, log=log)

        if merge_config:
            update_config(in_config=configfile, merge=merge_config, log=log)

        snake_command += ["--configfile", configfile]

        # display the runtime configuration
        snake_config = read_config(configfile)
        msg_box(
            "Runtime config",
            errmsg=yaml.dump(snake_config, Dumper=yaml.Dumper),
            log=log,
        )

    # add threads
    if "--profile" not in snake_args and profile is None:
        snake_command += ["--cores", threads]

    # add conda args if using conda
    if use_conda:
        snake_command += ["--use-conda"]
        if conda_prefix:
            snake_command += ["--conda-prefix", conda_prefix]
    if use_singularity:
        snake_command += ["--use-singularity"]
        if singularity_prefix:
            snake_command += ["--singularity-prefix", singularity_prefix]
    # add snakemake default args
    if snake_default:
        snake_command += snake_default

    # add any additional snakemake commands
    if snake_args:
        snake_command += list(snake_args)

    # allow double-handling of --profile
    if profile:
        snake_command += ["--profile", profile]

    # allow double-handling of --workflow-profile
    if workflow_profile:
        # copy system default if not present
        copy_config(
            os.path.join(workflow_profile, "config.yaml"),
            system_config=system_workflow_profile,
            log=log,
        )

        snake_command += ["--workflow-profile", workflow_profile]

    # Run Snakemake!!!
    snake_command = " ".join(str(s) for s in snake_command)
    msg_box("Snakemake command", errmsg=snake_command, log=log)
    if not subprocess.run(snake_command, shell=True).returncode == 0:
        msg("ERROR: Snakemake failed", log=log)
        sys.exit(1)
    else:
        msg("Snakemake finished successfully", log=log)
    return 0


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="zamp_out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads to use",
            default=1,
            show_default=True,
        ),
        click.option(
            "--use-singularity/--no-use-singularity",
            default=True,
            help="Use singularity containers for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--singularity-prefix",
            default=snake_base(os.path.join("singularity")),
            help="Custom singularity container directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=False,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("snakemake", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
            hidden=True,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.argument(
            "snake_args",
            nargs=-1,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def db_options(func):
    """
    Command line args for database command
    """
    options = [
        click.option(
            "--fasta",
            type=click.Path(readable=True),
            help="Path to database fasta file",
            required=True,
        ),
        click.option(
            "--taxonomy",
            type=click.Path(readable=True),
            help="Path to tab seperated taxonomy file in QIIME format",
            required=True,
        ),
        click.option(
            "--name",
            type=str,
            help="Database name",
            required=True,
        ),
        click.option(
            "--processing/--no-processing",
            default=True,
            help="Extract amplicon regions and merge taxonomy",
            show_default=True,
        ),
        click.option(
            "--tax-collapse",
            help="Dictionary of number of ranks to print limit when collapsing names ",
            default='{"Species": 5, "Genus": 6}',
            callback=validate_dict,
            show_default=True,
        ),
        click.option(
            "--fw-primer",
            type=str,
            help="Forward primer sequence to extract amplicon",
        ),
        click.option(
            "--rv-primer",
            type=str,
            help="Reverse primer sequence to extract amplicon",
        ),
        click.option(
            "--minlen",
            type=int,
            help="Minimum amplicon length",
            default=300,
            show_default=True,
        ),
        click.option(
            "--maxlen",
            type=int,
            help="Maximum amplicon length",
            default=500,
            show_default=True,
        ),
        click.option(
            "--ampcov",
            type=float,
            help="Minimum amplicon coverage",
            default=0.9,
            show_default=True,
        ),
        click.option(
            "--max-mismatch",
            type=int,
            help="Maximum number of accepted primer mismatches",
            default=4,
            show_default=True,
        ),
        click.option(
            "--rdp-mem",
            type=str,
            help="Maximum RAM for RDP training",
            default="30g",
            show_default=True,
        ),
        click.option(
            "--classifier",
            multiple=True,
            type=click.Choice(
                ["rdp", "qiimerdp", "dada2rdp", "decipher"], case_sensitive=False
            ),
            default=["rdp", "qiimerdp", "dada2rdp"],
            help="Which classifiers to train on the database",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func
