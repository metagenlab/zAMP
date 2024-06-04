### Match reads with sample list. By default will take the first column as Sample ID. If OldSampleName present, this column will be used instead to match the rows with the reads. They will be then renamed by the content of the first column.
import pandas
import re
import itertools


## Define a function to extract the .fastq name pattern
def get_read_naming_patterns(directory):
    result = []
    extension = {}
    for fname in sorted(os.listdir(directory)):
        if (
            fname.endswith("fastq.gz")
            or fname.endswith("fq.gz")
            or fname.endswith("fastq")
            or fname.endswith("fq")
        ):
            regex_str = r"(_L0+[1-9]+)?_(R)?(1|2)(\.|_)"  # regex for finding R1 and R2, if L001 is present before, it is also included
            regex = re.compile(regex_str)
            ext = re.search(regex, fname)
            if ext is None:
                ext = re.search(r"f(?:ast)?q(?:\.gz)?", fname)
                samp = re.sub(r'"\.$"', "", re.search(r"^([^\.]*)\.*", fname).group(0))
                if samp in extension.keys():
                    if ext.group(0).endswith(".gz"):
                        extension[samp] = [ext.group(0)]
                else:
                    extension[samp] = [ext.group(0)]
            else:
                regex_after = re.compile(regex_str + ".*")
                regex_before = re.compile(".*" + regex_str)
                read = re.compile(re.search(regex_after, fname).group(0))
                samp = re.sub(regex, "", re.search(regex_before, fname).group(0))
                extension.setdefault(samp, [])
                extension[samp].append(re.sub("^_", "", read.pattern))
    return extension


## Check for the presence of a directory for .fastq in config. If none, will be "links".
if "link_directory" in config.keys():
    link_directory = config["link_directory"]
    if not link_directory.endswith("/"):
        link_directory = link_directory + "/"
else:
    link_directory = "links/"

## Check that "tax_DB_path" ends with "/". Otherwise, correct.
if not config["tax_DB_path"].endswith("/"):
    config["tax_DB_path"] = config["tax_DB_path"] + "/"

## Check for the existence/accessbility of tax_DB_path
if os.path.isdir(config["tax_DB_path"]) is False:
    raise IOError("Cannot access to '{}' variable.".format("tax_DB_path"))


## Check for the presence of a metadata table in in config, either for local reads ("local_samples") or distant reads ("sra_samples")
if "local_samples" not in config.keys() and "sra_samples" not in config.keys():
    raise ValueError(
        "No samples defined in the config file. local_samples or sra_samples must be defined in config"
    )

## Check the local_samples/sra_samples tables for the presence of following variables in config: 'sample_label', 'grouping_column', 'run_column'.
### read sample dataframe from config file
if "local_samples" in config.keys():
    metadata = config["local_samples"]
elif "sra_samples" in config.keys():
    metadata = config["sra_samples"]
df = pandas.read_csv(metadata, sep="\t", index_col=0)
df_colnames = set(df.columns)

### Get values for which presence should be tested in dataframe
to_check = [
    "sample_label",
    "grouping_column",
    "run_column",
]  # a list of config values describing columns that should be in the dataframe
to_check_dict = dict((k, config[k]) for k in to_check)


### Set a function to compare list
def list_diff(list1, list2):
    return list(set(list1) - set(list2))


### Loop over the values of config variable to check
for key in to_check_dict.keys():

    ## Extract values from dictionnary
    to_check_values = to_check_dict[key]

    ## If it is a single value, test that it is in the columns of the metadata
    if isinstance(to_check_values, str):

        ## if the column name list is not in the column names of df dataframe then raise an error
        if to_check_values not in df_colnames:
            message = (
                f"{' and '.join(to_check_values)} not available in the sample TSV file."
            )
            raise IOError(message)

    ## If it is a list, check that it contained in metadata
    elif isinstance(to_check_values, list):

        diff = []
        diff = list_diff(to_check_values, df_colnames)

        ## if the set of column names in to_check list is not in the column names of df dataframe then raise an error
        if diff:
            for item in diff:
                """if the set of column names in to_check list is not a subset of column names of df dataframe then raise an error"""
                message = f"{' and '.join(diff)} not available in the sample TSV file."
                raise IOError(message)

    else:
        message = f"{key} should be a string or a list of string"
        raise IOError(message)


## Generate a list of input files, either local (when "local_samples" is in config, or SRA-hosted (when "sra_samples" is in config))
### Create the final pandas dataframe that will be completed
all_samples = pandas.DataFrame()

### Create a set of dictionaries to store sample characteristics
sras_ext = {}
reads_sra = {}
reads_local = {}
original_names = {}
reads_ext = {}
paths = {}
layout = {}


### In case of local samples, work our way through the local_samples table to extract read paths (if indicated in the R1/R2 columns) or extract match it with .fastq found in the "links" directory.
if "local_samples" in config.keys():
    ## Read the metadata table
    local_data = pandas.read_csv(config["local_samples"], sep="\t", index_col=0)
    local_data.index = [str(x) for x in local_data.index]
    all_local_sample_names = "".join(list(local_data.index))
    ## Check for forbidden characters in sample names
    if (
        "(" in all_local_sample_names
        or ")" in all_local_sample_names
        or "_-_" in all_local_sample_names
    ):
        raise ValueError("Forbidden character in sample name in sample name file")
    ## Extract path if indicated in "R1" column
    if "R1" in list(local_data):
        for sample, sample_data in local_data.iterrows():
            if sample in paths:
                ## Check for sample names to be unique
                raise IOError(f"Identical sample name used multiple times: {sample}")
            paths[sample] = [sample_data.loc["R1"]]
            reads_ext[sample] = "single"
            layout[sample] = "single"
            if "R2" in local_data.columns.values:
                if "R1" in str(sample_data.loc["R2"]):
                    raise IOError(
                        "ATTENTION! R1 flag within R2 filename: %s",
                        sample_data.loc["R2"],
                    )
                if (str(sample_data.loc["R2"]).strip()) != "nan":
                    paths[sample].append(sample_data.loc["R2"])
                    reads_ext[sample] = ["R1", "R2"]
                    layout[sample] = "paired"
        all_samples = local_data
        paths = {**paths}

    ## In absence of a "R1" column indicating file paths, try to match the sample names with the content of the "links" folder
    else:

        reads_local = get_read_naming_patterns(link_directory)
        original_names = {x: x for x in reads_local.keys()}
        read_correct = {}
        original_correct = {}

        ## In absence of a "OldSampleName", the leftmost column is used to match with "links" content
        if "OldSampleName" not in list(local_data):

            for sample in list(local_data.index):
                regex = re.compile(
                    r"%s([^a-zA-Z0-9]|$)" % sample
                )  # this regex ensures that the matching of the sample names end at the end of the str, to prevent S1 matching S10 for instance
                match = [
                    bool(re.match(regex, x))
                    for x in sorted(list(original_names.keys()))
                ]
                if (
                    sum(match) != 1
                ):  # there must be one and only one entry matching one sample name
                    raise ValueError(
                        "Problem matching SampleName to read file names: " + sample
                    )
                read_name = str(
                    sorted(list(original_names.keys()))[match.index(True)]
                )  # Sample name with _SX
                original_correct[sample] = original_names[read_name]
                read_correct[sample] = reads_local[read_name]
                paths[sample] = expand(
                    link_directory + read_name + "_{reads}", reads=read_correct[sample]
                )

                ## In presence of a "LibraryLayout" column, samples can be specifid to be "single" or "paired".
                if "LibraryLayout" in list(local_data):
                    if local_data.loc[sample, "LibraryLayout"].lower() == "paired":
                        reads_ext[sample] = ["R1", "R2"]
                        layout[sample] = "paired"
                    elif local_data.loc[sample, "LibraryLayout"].lower() == "single":
                        reads_ext[sample] = ["single"]
                        layout[sample] = "single"
                    else:
                        raise ValueError(
                            "Problem in the Local_sample file, LibraryLayout badly defined"
                        )
                ## Otherwise, "paired" is assumed
                else:
                    reads_ext[sample] = ["R1", "R2"]
                    layout[sample] = "paired"

        ## In presence of a "OldSampleName", the this column is used to match with "links" content
        else:
            for Old in list(local_data["OldSampleName"]):
                regex = re.compile(r"%s([^a-zA-Z0-9]|$)" % Old)
                match = [
                    bool(re.match(regex, x))
                    for x in sorted(list(original_names.keys()))
                ]
                if sum(match) != 1:
                    raise ValueError(
                        "Problem matching OldSampleName to read file names : " + Old
                    )
                read_name = str(
                    sorted(list(original_names.keys()))[match.index(True)]
                )  # Sample with have SX to unify with above
                sample = str(local_data.index[local_data["OldSampleName"] == Old][0])
                original_correct[sample] = original_names[read_name]
                read_correct[sample] = reads_local[read_name]
                paths[sample] = expand(
                    link_directory + read_name + "_{reads}", reads=read_correct[sample]
                )

                ## In presence of a "LibraryLayout" column, samples can be specifid to be "single" or "paired".
                if "LibraryLayout" in list(local_data):
                    if local_data.loc[sample, "LibraryLayout"].lower() == "paired":
                        reads_ext[sample] = ["R1", "R2"]
                        layout[sample] = "paired"
                    elif local_data.loc[sample, "LibraryLayout"].lower() == "single":
                        reads_ext[sample] = ["single"]
                        layout[sample] = "single"
                    else:
                        raise ValueError(
                            "Problem in the Local_sample file, LibraryLayout badly defined"
                        )
                ## Otherwise, "paired" is assumed
                else:
                    reads_ext[sample] = ["R1", "R2"]
                    layout[sample] = "paired"

        original_names = original_correct
        reads_local = read_correct
        reads_ext = reads_ext

        all_samples = local_data


### In case of sra-hosted samples, work our way through the "sra_table" to extract LibraryLayout and clean sample names.
if "sra_samples" in config.keys():
    sra_data = pandas.read_csv(
        config["sra_samples"], sep="\t", index_col=0
    ).drop_duplicates()
    all_sra_sample_names = "".join(list(sra_data.index))
    if (
        "(" in all_sra_sample_names
        or ")" in all_sra_sample_names
        or "_-_" in all_sra_sample_names
    ):
        raise ValueError("Forbidden character in sample name in sra file")
    for sra_sample in list(sra_data.index):
        sample_name = (
            str(sra_sample).replace(" ", "_").replace("&", "and").replace(":", "-")
        )
        if (
            sample_name in reads_sra.keys()
        ):  # if the sample name is already used, add _(n+1) at the end
            sample_name = (
                sample_name + "_" + str(list(reads_sra.keys()).count(sample_name))
            )
        reads_sra[sample_name] = str(sra_sample)
        if sra_data.loc[sra_sample, "LibraryLayout"].lower() == "paired":
            sras_ext[sample_name] = ["1.fastq.gz", "2.fastq.gz"]
            reads_ext[sample_name] = ["R1", "R2"]
            layout[sample_name] = "paired"
            paths[sample_name] = expand(
                link_directory + sample_name + "_{reads}", reads=sras_ext[sample_name]
            )
        elif sra_data.loc[sra_sample, "LibraryLayout"].lower() == "single":
            sras_ext[sample_name] = ["fastq.gz"]
            reads_ext[sample_name] = ["single"]
            layout[sample_name] = "single"
            paths[sample_name] = expand(
                link_directory + sample_name + ".{reads}", reads=sras_ext[sample_name]
            )
        else:
            raise ValueError("Problem in the sra file, LibraryLayout badly defined")
    all_samples = sra_data
    config["local_samples"] = config["sra_samples"]

    # all_samples.loc[sample_name, "Replicate"]=sra_data.loc[i, "Replicate"]
read_naming = {**reads_local, **sras_ext}
original_names = {**original_names, **reads_sra}
