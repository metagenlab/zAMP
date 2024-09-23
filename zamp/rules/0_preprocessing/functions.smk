# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


def get_raw_reads(wildcards):
    if SAMPLES.loc[wildcards.sample, "paired"]:
        return list(SAMPLES.loc[wildcards.sample, "R1":"R2"])
    else:
        return list(SAMPLES.loc[wildcards.sample, "R1"])
