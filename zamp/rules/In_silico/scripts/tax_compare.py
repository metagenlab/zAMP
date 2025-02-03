import pandas as pd
import os
import numpy as np
from snakemake.script import snakemake

# Variables
ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
# Make plaheloder name dict
placeholder = {}
for rank in ranks:
    placeholder[f"{rank}"] = "_" + rank[0].lower()


# Functions
def get_basename(path, suffix):
    return os.path.basename(path).split(suffix)[0].rstrip(".")


def propagate_nan(df):
    nan_mask = df.isna()
    nan_cumsum = nan_mask.cumsum(axis=1)
    df[nan_cumsum > 0] = np.nan
    return df


def get_full_species_names(species):
    if "(" in species:
        species = species.split("(")[0]
    if "/" in species:
        species_list = species.split("/")
        newnames = []
        for sp in species_list:
            if (" " not in sp) and ("_s" not in sp):
                sp = species_list[0].split(" ")[0] + " " + sp
                newnames.append(sp)
            else:
                newnames.append(sp)
        return "/".join(newnames)
    else:
        return species


def get_match(expected, assigned, rank):
    if "/" in assigned:
        return expected in assigned
    else:
        if f"_{rank[0]}" in assigned:
            return False
        else:
            return expected == assigned


## Read count table
count_df = (
    pd.read_csv(snakemake.input.count_table, sep="\t")
    .transpose()
    .melt(var_name="seq_id", value_name="seq_count", ignore_index=False)
)
count_df.index.name = "fasta"
count_df.reset_index(inplace=True)

## Read genome assembly taxonomy
if snakemake.params.local:
    asm_df = pd.read_csv(snakemake.input.expected_taxonomy, sep="\t")
    asm_df.fasta = asm_df.fasta.apply(
        lambda x: get_basename(x, snakemake.params.suffix)
    )
else:
    asm_df = pd.read_csv(snakemake.input.expected_taxonomy[0], sep="\t")
    asm_df.fasta = asm_df.path.apply(lambda x: get_basename(x, snakemake.params.suffix))
    # asm_df.fasta = asm_df.accession + "_" + asm_df.fasta
    asm_df = asm_df[["accession", "fasta", "assembly_level"]].merge(
        pd.read_csv(snakemake.input.expected_taxonomy[1], sep="\t"), on="accession"
    )

## Read database taxonomy
db_tax_df = pd.read_csv(snakemake.input.db_tax, sep="\t", names=["seq_id", "tax"])
db_tax_df[ranks] = db_tax_df.tax.str.split(";", expand=True).loc[:, 0:6]
db_tax_df["full_species_names"] = db_tax_df.species.apply(get_full_species_names)
asm_df["genus_in_DB"] = asm_df.genus.isin(db_tax_df.genus)
asm_df["species_in_DB"] = asm_df["species"].apply(
    lambda x: any(x in match for match in db_tax_df["full_species_names"])
)

## Read assigned taxonomy
assigned_tax_df = pd.read_csv(snakemake.input.assigned_taxonomy, sep="\t", header=None)
assigned_tax_df.columns = ["seq_id", "tax", "confidence"]
assigned_tax_df[ranks] = assigned_tax_df.tax.str.split(";", expand=True)

## Add taxonomy to sequences in count table
amp_df = (
    count_df[count_df.seq_count > 0]
    .merge(assigned_tax_df, on="seq_id", how="left")
    .sort_values("fasta")
)
## Replace no amp sequence counts by 0
amp_df.loc[amp_df["seq_id"] == "No_amp", "seq_count"] = 0
amp_df["amplified"] = amp_df.seq_id.apply(lambda x: False if x == "No_amp" else True)
amp_df = amp_df.merge(asm_df, on="fasta", suffixes=("_assigned", "_expected"))

amp_df["assigned_species_full"] = amp_df.species_assigned.apply(
    lambda x: get_full_species_names(x) if pd.notna(x) else x
)
amp_df["genus_match"] = amp_df.apply(
    lambda row: get_match(
        str(row["genus_expected"]), str(row["genus_assigned"]), rank="genus"
    ),
    axis=1,
)
amp_df["species_match"] = amp_df.apply(
    lambda row: get_match(
        str(row["species_expected"]), str(row["assigned_species_full"]), rank="species"
    ),
    axis=1,
)
# Add summary stats in assembly_table
asm_df = asm_df.merge(
    amp_df.groupby(["fasta", "amplified"], as_index=False).agg(
        {"seq_count": "sum", "genus_match": "mean", "species_match": "mean"}
    ),
    on="fasta",
)
# asm_df["expected_tax"] = asm_df[ranks].apply(lambda x: ";".join(x), axis=1)
# asm_df.drop(ranks, axis=1, inplace=True)

# Save tables
amp_df.to_csv(snakemake.output.amplicons, sep="\t", index=False)
asm_df.to_csv(snakemake.output.assemblies, sep="\t", index=False)
