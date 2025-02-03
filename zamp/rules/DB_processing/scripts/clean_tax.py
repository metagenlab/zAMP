import pandas as pd
import numpy as np
from snakemake.script import snakemake


def propagate_nan(df):
    # Identify columns that have NaN values
    nan_mask = df.isna()

    # Create a cumulative sum mask to propagate NaNs
    nan_cumsum = nan_mask.cumsum(axis=1)

    # Any value where the cumulative sum is greater than 0 should be NaN
    df[nan_cumsum > 0] = np.nan

    return df


df = pd.read_csv(snakemake.input[0], sep="\t")
if any(df.columns.str.contains(";")):
    df = pd.read_csv(snakemake.input[0], sep="\t", header=None)

df.columns = ["seq_id", "tax"]

ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# Make plaheloder name dict
placeholder = {}
for rank in ranks:
    placeholder[f"{rank}"] = "_placeholder_" + rank[0].lower()

if "unite" in snakemake.params.db_name:
    df.tax = df.tax.replace(to_replace=r"[a-z]__", value="", regex=True)
    df.tax = df.tax.replace(to_replace=r"_", value=" ", regex=True)

if "greengenes" in snakemake.params.db_name:
    df.tax = df.tax.str.replace("; ", ";")
    ## Remove leading k__ to s__ in GTDB taxonomy
    df.tax = df.tax.replace(to_replace=r"[a-z]__", value="", regex=True)

lintax_df = df.tax.str.split(";", expand=True).loc[:, 0:6]
lintax_df.columns = ranks
lintax_df = lintax_df.replace("", np.nan).fillna(np.nan)

if "silva" in snakemake.params.db_name:
    ## Replace taxa containing and Unkown or Incertae with NaN
    lintax_df.replace(
        to_replace=r".*Incertae.*|.*Unknown.*", value=np.nan, regex=True, inplace=True
    )

    ## Replace endosymbionts by NaN
    lintax_df.replace("endosymbionts", np.nan, inplace=True)


if "greengenes" in snakemake.params.db_name:
    # In greeengenes2, some sequences have the same taxa name assigned to multiple ranks
    # the code below iterates over ranks and replaces duplicate names with NaN
    # This mitigates convergent taxonomy errors in RDP
    array = lintax_df.to_numpy()

    # Iterate through the array to replace duplicates with NaN
    for i in range(array.shape[0]):
        seen = set()
        for j in range(array.shape[1]):
            if array[i, j] in seen:
                array[i, j] = np.nan
            else:
                seen.add(array[i, j])
    lintax_df = pd.DataFrame(array, columns=lintax_df.columns)


filled_df = propagate_nan(lintax_df).ffill(axis=1)
prop_tax_df = pd.DataFrame()
for n, rank in enumerate(ranks):
    if rank != ranks[0]:
        prev = ranks[n - 1]
        duplicated = (
            filled_df[filled_df[f"{prev}"] == filled_df[f"{rank}"]][f"{rank}"]
            + f"{placeholder[rank]}"
        )
        prop_tax_df[f"{rank}"] = lintax_df[f"{rank}"].combine_first(duplicated)
    else:
        prop_tax_df[f"{rank}"] = filled_df[f"{rank}"]

if "silva" in snakemake.params.db_name:
    ## Get classified species index
    index = ~prop_tax_df["Species"].str.contains(placeholder["Species"])
    ## Add genus name in species for classified species
    prop_tax_df.loc[index, "Species"] = (
        prop_tax_df.loc[index, "Genus"] + " " + prop_tax_df.loc[index, "Species"]
    )

prop_tax_df["taxpath"] = prop_tax_df[ranks].T.agg(";".join)
prop_tax_df.taxpath = prop_tax_df.taxpath.str.replace("_placeholder", "")
df = df.join(prop_tax_df)
df[["seq_id", "taxpath"]].to_csv(
    snakemake.output[0], sep="\t", index=False, header=False
)
