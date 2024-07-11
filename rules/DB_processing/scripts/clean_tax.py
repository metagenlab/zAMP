import pandas as pd
import numpy as np


def propagate_nan(df):
    # Identify columns that have NaN values
    nan_mask = df.isna()

    # Create a cumulative sum mask to propagate NaNs
    nan_cumsum = nan_mask.cumsum(axis=1)

    # Any value where the cumulative sum is greater than 0 should be NaN
    df[nan_cumsum > 0] = np.nan

    return df


df = pd.read_csv(snakemake.input[0], sep="\t")
df.columns = ["seq_id", "tax"]

ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

if snakemake.params.db_version == "greengenes2":
    df.tax = df.tax.str.replace("; ", ";")
    ## Remove leading k__ to s__ in GTDB taxonomy
    df.tax = df.tax.replace(to_replace=r"[a-z]__", value="", regex=True)

lintax_df = df.tax.str.split(";", expand=True).loc[:, 0:6]
lintax_df.columns = ranks
lintax_df = lintax_df.replace("", np.nan).fillna(np.nan)


if snakemake.params.db_version == "greengenes2":
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
        if rank == "Species":
            placeholder = rank[0:2].lower()
        else:
            placeholder = rank[0].lower()
        duplicated = (
            filled_df[filled_df[f"{prev}"] == filled_df[f"{rank}"]][f"{rank}"]
            + f"_{placeholder}"
        )
        prop_tax_df[f"{rank}"] = lintax_df[f"{rank}"].combine_first(duplicated)
    else:
        prop_tax_df[f"{rank}"] = filled_df[f"{rank}"]

prop_tax_df["taxpath"] = prop_tax_df[ranks].T.agg(";".join)
df = df.join(prop_tax_df)
df[["seq_id", "taxpath"]].to_csv(
    snakemake.output[0], sep="\t", index=False, header=False
)
