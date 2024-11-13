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
def get_basename(path):
    return os.path.basename(path).split("_genomic.fna")[0]


def propagate_nan(df):
    nan_mask = df.isna()
    nan_cumsum = nan_mask.cumsum(axis=1)
    df[nan_cumsum > 0] = np.nan
    return df


## Read count table
count_df = (
    pd.read_csv(snakemake.input.count_table, sep="\t")
    .transpose()
    .melt(var_name="seq_id", value_name="seq_count", ignore_index=False)
)
count_df.index.name = "assembly_name"
count_df.reset_index(inplace=True)

## Read genome assembly metadata with their taxonomy
asm_df = pd.read_csv(snakemake.input.assembly_summary, sep="\t")
asm_df.assembly_name = asm_df.path.apply(get_basename)
asm_df.assembly_name = asm_df.accession + "_" + asm_df.assembly_name
asm_df = asm_df[["accession", "assembly_name", "assembly_level"]].merge(
    pd.read_csv(snakemake.input.expected_taxonomy, sep="\t"), on="accession"
)
## Read database taxonomy
db_tax_df = pd.read_csv(snakemake.input.db_tax, sep="\t")
asm_df["in_DB"] = asm_df.species.isin(db_tax_df.Species)
asm_df.set_index("assembly_name", inplace=True)

## Read assigned taxonomy
assigned_tax_df = pd.read_csv(snakemake.input.assigned_taxonomy, sep="\t", header=False)
assigned_tax_df.columns = ["seq_id", "tax", "confidence"]
tax = assigned_tax_df.tax.str.split(";", expand=True)
tax.columns = ranks
### Add placeholder names for unassigned ranks
filled_tax = propagate_nan(tax).ffill(axis=1)
fixed_tax = pd.DataFrame()
for n, rank in enumerate(ranks):
    if rank != ranks[0]:
        prev = ranks[n - 1]
        duplicated = (
            filled_tax[filled_tax[f"{prev}"] == filled_tax[f"{rank}"]][f"{rank}"]
            + f"{placeholder[rank]}"
        )
        fixed_tax[f"{rank}"] = tax[f"{rank}"].combine_first(duplicated)
    else:
        fixed_tax[f"{rank}"] = filled_tax[f"{rank}"]
assigned_tax_df = assigned_tax_df.join(fixed_tax)

## Add taxonomy to sequences in count table
amp_df = (
    count_df[count_df.seq_count > 0]
    .merge(assigned_tax_df, on="seq_id", how="left")
    .sort_values("assembly_name")
)
## Replace no amp sequence counts by 0
amp_df.loc[amp_df["seq_id"] == "No_amp", "seq_count"] = 0
amp_df["genus_match"] = amp_df.genus.isin(asm_df.genus)
amp_df["species_match"] = amp_df.species.isin(asm_df.species)
amp_df["amplified"] = amp_df.seq_id.apply(lambda x: False if x == "No_amp" else True)

# Add summary stats in assembly_table
asm_df = asm_df.merge(
    amp_df.groupby(["assembly_name", "amplified"], as_index=False).agg(
        {"seq_count": "sum", "genus_match": "mean", "species_match": "mean"}
    ),
    on="assembly_name",
)
asm_df.drop(ranks[:-2], axis=1, inplace=True)

# Save tables
amp_df.to_csv(snakemake.output.amplicons, sep="\t", index=False)
asm_df.to_csv(snakemake.output.assemblies, sep="\t", index=False)
