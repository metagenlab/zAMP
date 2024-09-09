# dirs
dir = AttrMap()

# Base output dir
try:
    assert (au.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "zamp_out"

# Output subdirs
## Reads
dir.out.fastq = os.path.join(dir.out.base, "fastq")
## QC
dir.out.qc = os.path.join(dir.out.base, "QC")
## Trimming
dir.out.cutadapt = os.path.join(dir.out.base, "cutadapt")
dir.out.pandaseq = os.path.join(dir.out.base, "pandaseq")
## Denoising
dir.out.vsearch = os.path.join(dir.out.base, "vsearch")
dir.out.dada2 = os.path.join(dir.out.base, "DADA2")
# logs
dir.logs = os.path.join(dir.out.base, "logs")
# envs
dir.envs = os.path.join(workflow.basedir, "envs")
