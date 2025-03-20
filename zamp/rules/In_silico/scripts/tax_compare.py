import pandas as pd
import os
from snakemake.script import snakemake

# Variables
ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


# Functions
def get_basename(path, suffix):
    return os.path.basename(path).split(suffix)[0].rstrip(".")


def get_match(expected, assigned):
    return all(word in str(assigned) for word in str(expected).split())


def get_inconsitent_duplicates(df):
    duplicates = df[df.duplicated("formatted_tax", keep=False)]
    return duplicates.groupby("formatted_tax").filter(
        lambda x: x["all_tax"].nunique() > 1
    )


## Read count table
count_df = pd.read_csv(snakemake.input.count_table, sep="\t")


## Read genome assembly taxonomy
if snakemake.params.local:
    asm_df = pd.read_csv(snakemake.input.expected_taxonomy, sep="\t")
    asm_df.fasta = asm_df.fasta.apply(
        lambda x: get_basename(x, snakemake.params.suffix)
    )
else:
    asm_df = pd.read_csv(snakemake.input.expected_taxonomy[0], sep="\t")
    asm_df.fasta = asm_df.path.apply(lambda x: get_basename(x, snakemake.params.suffix))
    asm_df = asm_df[["accession", "fasta", "assembly_level"]].merge(
        pd.read_csv(snakemake.input.expected_taxonomy[1], sep="\t"), on="accession"
    )

## Read database taxonomy
db_tax_df = pd.read_csv(snakemake.input.db_tax, sep="\t")
db_tax_df[ranks] = db_tax_df.all_tax.str.split(";", expand=True).loc[:, 0:6]
asm_df["genus_in_DB"] = asm_df["genus"].apply(
    lambda x: any(x in match for match in db_tax_df["genus"])
)
asm_df["species_in_DB"] = asm_df["species"].apply(
    lambda x: any(x in match for match in db_tax_df["species"])
)

## Read assigned taxonomy
assigned_tax_df = pd.read_csv(snakemake.input.assigned_taxonomy, sep="\t", header=None)
if "formatted_tax" in db_tax_df.columns:
    assigned_tax_df.columns = ["asv_id", "formatted_tax", "confidence"]
    # Get formatted to full taxonomy dictionary
    tax2all = dict(zip(db_tax_df.formatted_tax, db_tax_df.all_tax))
    # Get all_tax column
    all_tax = []
    for tax in assigned_tax_df.formatted_tax:
        try:
            all_tax.append(tax2all[tax])
        except KeyError:
            all_tax.append(None)
    assigned_tax_df["all_tax"] = all_tax
    assigned_tax_df["all_tax"] = assigned_tax_df["all_tax"].fillna(
        assigned_tax_df["formatted_tax"]
    )

else:
    assigned_tax_df.columns = ["asv_id", "all_tax", "confidence"]

assigned_tax_df[ranks] = assigned_tax_df.all_tax.str.split(";", expand=True)

## Add taxonomy to sequences in count table
amp_df = count_df.merge(assigned_tax_df, on="asv_id", how="left").sort_values("fasta")
amp_df["amplified"] = amp_df.asv_id.apply(lambda x: False if x == "no_amp" else True)
amp_df = amp_df.merge(asm_df, on="fasta", suffixes=("_assigned", "_expected"))


amp_df["genus_match"] = amp_df.apply(
    lambda row: get_match(
        row["genus_expected"],
        row["genus_assigned"],
    ),
    axis=1,
)
amp_df["species_match"] = amp_df.apply(
    lambda row: get_match(
        row["species_expected"],
        row["species_assigned"],
    ),
    axis=1,
)
# Add summary stats in assembly_table
asm_df = asm_df.merge(
    amp_df.groupby(["fasta", "amplified"], as_index=False).agg(
        {"asv_count": "sum", "genus_match": "mean", "species_match": "mean"}
    ),
    on="fasta",
)

# Save tables
amp_df.to_csv(snakemake.output.amplicons, sep="\t", index=False)
asm_df.to_csv(snakemake.output.assemblies, sep="\t", index=False)
