import os
import click
import yaml
import sys
import subprocess


from snaketool_utils.cli_utils import (
    msg,
    read_config,
    update_config,
    copy_config,
    echo_click,
    msg_box,
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
    use_singularity=True,
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

    # add snakemake default args
    if snake_default:
        snake_command += snake_default

    # add any additional snakemake commands
    if snake_args:
        snake_command += list(snake_args)

    # allow double-handling of --profile
    if profile:
        snake_command += ["--profile", profile]
    if use_singularity:
        snake_command += ["--use-singularity"]
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
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
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
            help="Comma seperated list of database names",
            required=True,
        ),
        click.option(
            "--processing/--no-processing",
            default=True,
            help="Extract amplicon regions and merge taxonomy",
            show_default=True,
        ),
        click.option(
            "--fw-primer",
            type=str,
            help="Forward primer sequence to extract amplicon",
            required=True,
        ),
        click.option(
            "--rv-primer",
            type=str,
            help="Reverse primer sequence to extract amplicon",
            required=True,
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
            "--errors",
            help="Maximum number of accepted primer mismatches, or float between 0 and 1",
            default=0.1,
            show_default=True,
        ),
        click.option(
            "--cutadapt_args_fw",
            type=str,
            default="",
            help="Additional cutadapt arguments for forward primer",
        ),
        click.option(
            "--cutadapt_args_rv",
            type=str,
            default="",
            help="Additional cutadapt arguments for reverse primer",
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
                ["rdp", "qiime2", "dada2", "decipher", "sintax"], case_sensitive=False
            ),
            default=["rdp", "qiime2", "dada2", "sintax"],
            help="Which classifiers to train on the database",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def run_options(func):
    """
    Command line args for running zAMP
    """
    options = [
        click.option(
            "--input",
            "-i",
            type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
            help="Input file/directory for fastq path",
            required=True,
        ),
        click.option(
            "--local/--sra",
            default=True,
            help="Whether reads are local or to be downloaded from NCBI SRA",
            show_default=True,
        ),
        click.option(
            "--metadata",
            "-m",
            type=click.Path(exists=True, file_okay=True, readable=True),
            help="Path to tab seperated metadata file",
        ),
        click.option(
            "--database",
            "-db",
            type=click.Path(exists=True, dir_okay=True, readable=True),
            help="Path to Database directory",
            required=True,
        ),
        click.option(
            "--name",
            type=str,
            help="Comma seperated list of database names",
        ),
        click.option(
            "--denoiser",
            multiple=True,
            type=click.Choice(["DADA2", "vsearch"], case_sensitive=False),
            default=["DADA2"],
            help="Choose dada2 or vsearch for denoising reads",
            show_default=True,
        ),
        click.option(
            "--classifier",
            multiple=True,
            type=click.Choice(
                ["rdp", "qiime2", "dada2", "decipher", "sintax"], case_sensitive=False
            ),
            default=["rdp"],
            help="Which classifiers to use for taxonomic assignment",
            show_default=True,
        ),
        click.option(
            "--trim/--no-trim",
            default=True,
            help="Trim primers or not",
            show_default=True,
        ),
        click.option(
            "--min-overlap",
            default=10,
            help="Minimum R1 and R2 overlap for reads merging",
            show_default=True,
        ),
        click.option(
            "--fw-primer",
            type=str,
            help="Forward primer sequence to extract amplicon from reads",
            required=True,
        ),
        click.option(
            "--rv-primer",
            type=str,
            help="Reverse primer sequence to extract amplicon from reads",
            required=True,
        ),
        click.option(
            "--minlen",
            default=390,
            type=int,
            help="Minimum read length for merged reads",
            show_default=True,
        ),
        click.option(
            "--maxlen",
            default=480,
            type=int,
            help="Maximum read length for merged reads",
            show_default=True,
        ),
        click.option(
            "--fw-trim",
            default=280,
            type=int,
            help="Minimum read length to trim low quality ends of R1 for DADA2 denoising",
            show_default=True,
        ),
        click.option(
            "--rv-trim",
            default=255,
            type=int,
            help="Minimum read length to trim low quality ends of R2 for DADA2 denoising",
            show_default=True,
        ),
        click.option(
            "--fw-errors",
            default=10,
            type=int,
            help="Maximum expected errors in R1 for DADA2 denoising",
            show_default=True,
        ),
        click.option(
            "--rv-errors",
            default=10,
            type=int,
            help="Maximum expected errors in R2 for DADA2 denoising",
            show_default=True,
        ),
        click.option(
            "--rarefaction",
            type=str,
            default="50000",
            help="Comma seperated list of number of reads for rarefaction",
            show_default=True,
        ),
        click.option(
            "--min-prev",
            default=0,
            type=float,
            help="Proporition (in %) of samples in which the feature has to be found to be kept",
            show_default=True,
        ),
        click.option(
            "--min-count",
            default=0,
            type=int,
            help="Minimal reads to be kept",
            show_default=True,
        ),
        click.option(
            "--normalization",
            type=str,
            default="NONE",
            help="Comma seperated list of values for counts normalization",
            show_default=True,
        ),
        click.option(
            "--replace-empty/--no-replace-empty",
            default=False,
            help="Replace empty taxa by placeholders",
            show_default=True,
        ),
        click.option(
            "--keep-rank",
            type=click.Choice(
                ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
            ),
            default=["Kingdom"],
            multiple=True,
            help="Rank to keep taxon",
            show_default=True,
        ),
        click.option(
            "--keep-taxa",
            type=str,
            default="Bacteria",
            help="Comma seperated list of taxa to keep",
            show_default=True,
        ),
        click.option(
            "--exclude-rank",
            type=str,
            default="Phylum",
            help="Comma seperated list of ranks to exclude",
            show_default=True,
        ),
        click.option(
            "--exclude-taxa",
            type=str,
            default="Bacteria_phy",
            help="Comma seperated list of taxa to exclude",
            show_default=True,
        ),
        click.option(
            "--melted/--no-melted",
            default=False,
            help="Generate melted phyloseq table",
            show_default=True,
        ),
        click.option(
            "--physeq-rank",
            type=str,
            help="Comma seperated list of ranks to collapse on in phyloseq output",
            default="OTU",
            show_default=True,
        ),
        click.option(
            "--transposed/--no-transposed",
            default=False,
            help="Transposed count table",
            show_default=True,
        ),
        click.option(
            "--qiime-viz",
            default=True,
            help="Output QIIME visualisation",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def insilico_options(func):
    """
    Command line args for database command
    """
    options = [
        click.option(
            "--input",
            "-i",
            type=str,
            help="Input fasta directory or taxa/accessions",
            required=True,
        ),
        click.option(
            "--tax",
            type=click.Path(file_okay=True, exists=True),
            help="Input taxonomic table when using local fasta",
        ),
        click.option(
            "--taxon/--accession",
            help="Are queries taxa or accessions",
            type=bool,
            default=True,
            show_default=True,
        ),
        click.option(
            "--limit",
            "-nb",
            help="Limit number of genomes per query",
            type=str,
            default=1,
        ),
        click.option(
            "--assembly",
            help="Comma seperated list of assembly level: complete,chromosome,scaffold,contig",
            default=None,
            show_default=True,
        ),
        click.option(
            "--only-ref/--not-only-ref",
            type=bool,
            help="Limit to reference and representative genomes",
            default=True,
            show_default=True,
        ),
        click.option(
            "--rank",
            help="taxonomic rank to filter by assemblies ",
            default=None,
            type=click.Choice(
                [
                    "superkingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                ],
                case_sensitive=False,
            ),
            show_default=True,
        ),
        click.option(
            "--nrank",
            help="Number of genomes per taxonomic rank",
            type=int,
            default=None,
            show_default=True,
        ),
        click.option(
            "--database",
            "-db",
            type=click.Path(exists=True, dir_okay=True, readable=True),
            help="Path to Database directory",
            required=True,
        ),
        click.option(
            "--name",
            type=str,
            help="Comma seperated list of database names",
        ),
        click.option(
            "--pcr-tool",
            type=click.Choice(["simulate", "in-silico"], case_sensitive=False),
            help="Tool for in silico PCR",
            default="in-silico",
            show_default=True,
        ),
        click.option(
            "--af-args",
            type=str,
            default="",
            help="Assembly_finder arguments",
        ),
        click.option(
            "--mismatch",
            help="Number of mismatches for simulate_PCR",
            type=int,
            default=3,
            show_default=True,
        ),
        click.option(
            "--threeprime",
            help="Number of match at the 3' end for a hit to be considered for simulate_PCR",
            type=int,
            default=2,
            show_default=True,
        ),
        click.option(
            "--fw-primer",
            type=str,
            help="Forward primer sequence to extract amplicon from reads",
            required=True,
        ),
        click.option(
            "--rv-primer",
            type=str,
            help="Reverse primer sequence to extract amplicon from reads",
            required=True,
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
            "--errors",
            help="Maximum number of accepted primer mismatches, or float between 0 and 1",
            default=0.1,
            show_default=True,
        ),
        click.option(
            "--classifier",
            type=str,
            default="rdp",
            help="classifier choise for taxonomic assignment [rdp, qiime2, dada2, sintax]",
            show_default=True,
        ),
        click.option(
            "--denoiser",
            multiple=True,
            type=click.Choice(["DADA2", "vsearch"], case_sensitive=False),
            default=["DADA2"],
            help="Choose dada2 or vsearch for denoising reads",
            show_default=True,
        ),
        click.option(
            "--replace-empty/--no-replace-empty",
            default=False,
            help="Replace empty taxa by placeholders",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func
