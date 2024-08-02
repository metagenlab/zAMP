def get_raw_reads(wildcards):
    if SAMPLES.loc[wildcards.sample, "paired"]:
        return list(SAMPLES.loc[wildcards.sample, "R1":"R2"])
    else:
        return list(SAMPLES.loc[wildcards.sample, "R1"])
