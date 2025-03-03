import pandas as pd
import numpy as np
from snakemake.script import snakemake


def find_convergent_taxa(df):
    """
    Identifies rows where a taxon is duplicated but has a different origin
    (Same species but two different genera for example)
    """
    inconsistent_rows = []

    # Iterate over columns starting from the second one
    for i in range(1, len(df.columns)):
        col = df.columns[i]
        prev_col = df.columns[i - 1]

        # Group by the current column to find duplicates
        grouped = df.groupby(col)

        for value, group in grouped:
            # If there are duplicates in the current column
            if len(group) > 1:
                # Check if the previous column values are different
                if group[prev_col].nunique() > 1:
                    inconsistent_rows.append(group)

    # Combine all inconsistent groups into a single DataFrame
    if inconsistent_rows:
        return pd.concat(inconsistent_rows).drop_duplicates()
    else:
        return pd.DataFrame()  # Return empty DataFrame if no inconsistencies are found


df = pd.read_csv(snakemake.input[0], sep="\t")
if any(df.columns.str.contains(";")):
    df = pd.read_csv(snakemake.input[0], sep="\t", header=None)

df.columns = ["seq_id", "tax"]

ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


if "unite" in snakemake.params.db_name:
    df.tax = df.tax.replace(to_replace=r"[a-z]__", value="", regex=True)
    df.tax = df.tax.replace(to_replace=r"_", value=" ", regex=True)

if "greengenes" in snakemake.params.db_name:
    df.tax = df.tax.str.replace("; ", ";")
    ## Remove leading k__ to s__ in GTDB taxonomy
    df.tax = df.tax.replace(to_replace=r"[a-z]__", value="", regex=True)

df[ranks] = df.tax.str.split(";", expand=True).loc[:, 0:6]


if "silva" in snakemake.params.db_name:
    ## Replace taxa containing and Unkown or Incertae with NaN
    df.replace(
        to_replace=r".*Incertae.*|.*Unknown.*", value=np.nan, regex=True, inplace=True
    )

    ## Replace endosymbionts by NaN
    df.replace("endosymbionts", np.nan, inplace=True)


if "greengenes" in snakemake.params.db_name:
    # In greeengenes2, some sequences have the same taxa name assigned to multiple ranks
    # the code below iterates over ranks and replaces duplicate names with NaN
    # This mitigates convergent taxonomy errors in RDP
    array = df[ranks].to_numpy()

    # Iterate through the array to replace duplicates with NaN
    for i in range(array.shape[0]):
        seen = set()
        for j in range(array.shape[1]):
            if array[i, j] in seen:
                array[i, j] = np.nan
            else:
                seen.add(array[i, j])
    df[ranks] = pd.DataFrame(array, columns=ranks)

df = df.replace("", np.nan).fillna(np.nan)
## Propagate ranks
nans = df.isna()
for n, rank in enumerate(ranks):
    if rank != "Kingdom":
        prev = ranks[n - 1]
        df.loc[nans[rank], rank] = (
            df.loc[nans[rank], prev] + "_placeholder_" + rank[0].lower()
        )

if "silva" in snakemake.params.db_name:
    ## Get classified species index
    index = ~df["Species"].str.contains("_placeholder_", na=False)
    ## Add parent name in species for classified species
    df.loc[index, "Species"] = df.loc[index, "Genus"] + " " + df.loc[index, "Species"]

conv_df = find_convergent_taxa(df[ranks])
if not conv_df.empty:
    raise ValueError(f"These rows cause convergent evolution:\n{print(conv_df)}")

df["all_tax"] = df[ranks].T.agg(";".join)
df.all_tax = df.all_tax.str.replace("_placeholder", "")

df[["seq_id", "all_tax"]].to_csv(
    snakemake.output[0], sep="\t", index=False, header=False
)
