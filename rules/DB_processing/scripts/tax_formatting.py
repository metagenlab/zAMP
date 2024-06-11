import pandas as pd
import numpy as np


# Functions
def dedup_set(s):
    return list(set(s))[0]


def format_discrepant_tax(rank, tax, rank_lim=None, print_desc=True):
    total_nb = len(tax)
    if rank_lim:
        try:
            nb_print = rank_lim[f"{rank}"]
            tax = list(tax)[0:nb_print]
        except KeyError:
            tax = list(tax)

    if total_nb > 1:
        if print_desc:
            prefix = f"Disc.{rank}_"
        else:
            prefix = ""

        return prefix + "/".join(tax) + f"({total_nb})"
    else:
        return list(tax)[0]


def problematic_taxa(row):
    return any("/" in str(cell) for cell in row)


# Ranks list
ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
ranks_lim = {
    "Species": snakemake.params.numbers_species,
    "Genus": snakemake.params.numbers_genus,
}

# Read tables
tax_df = pd.read_csv(snakemake.input.tax, sep="\t", names=["seq_id", "tax"])
uc_df = pd.read_csv(
    snakemake.input.uc,
    sep="\t",
    names=[
        "record_type",
        "cluster_id",
        "seq_length",
        "percent_identity",
        "strand",
        "not_used1",
        "not_used2",
        "comp_align",
        "query_id",
        "target_id",
    ],
)
# Taxonomy table
## Split taxonomy path into one column per rank
lintax_df = tax_df.tax.str.split(";", expand=True).loc[:, 0:6]
lintax_df.columns = ranks
## Replace empty values by NaN
lintax_df = lintax_df.replace("", np.nan).fillna(np.nan)

## Replace NaN values by latest no NaN value (forward fill)
filled_df = lintax_df.ffill(axis=1)

## Propagate latest non NaN taxon to empty ones
## Check repeated values and append them with first letter rank suffix
prop_tax_df = pd.DataFrame()
for n, rank in enumerate(ranks):
    if rank != "kingdom":
        prev = ranks[n - 1]
        duplicated = (
            filled_df[filled_df[f"{prev}"] == filled_df[f"{rank}"]][f"{rank}"]
            + f"_{rank[0].lower()}"
        )
        prop_tax_df[f"{rank}"] = lintax_df[f"{rank}"].combine_first(duplicated)
    else:
        prop_tax_df[f"{rank}"] = filled_df[f"{rank}"]

## Get classified species index
index = ~prop_tax_df["Species"].str.contains("_s")
## Add genus name in species for classified species
prop_tax_df.loc[index, "Species"] = (
    prop_tax_df.loc[index, "Genus"] + " " + prop_tax_df.loc[index, "Species"]
)


## Add seq_id to formatted taxonomy
tax_df = prop_tax_df.join(tax_df)

# Vsearch output table
## Extract C and H hits, check https://www.drive5.com/usearch/manual6/ucout.html for more info
clust_df = uc_df[uc_df.record_type.isin(["C", "H"])][
    ["record_type", "cluster_id", "query_id", "target_id"]
]

## Rename columns
clust_df.columns = ["record_type", "clust_id", "seq_id", "clust_rep"]
## Replace * with seqid
clust_df.loc[clust_df["clust_rep"] == "*", "clust_rep"] = clust_df["seq_id"]

## Add taxonomy to clusters
clust_df = clust_df.merge(tax_df, on="seq_id")

## Get unique taxonomy indexes
ranks_idxs = (
    clust_df.groupby(ranks)["clust_id"].transform(lambda x: pd.factorize(x)[0]) + 1
)
clust_df.insert(clust_df.shape[1], "rank_idx", ranks_idxs)


## For unclassified taxa, add index in taxon name for each cluster
### Example : Two Pseudoclavibacter_s sequences clusterin into clusters
###           will be named Pseudoclavibacter_s1 and Pseudoclavibacter_s2 for each cluster
for rank in ranks:
    clust_df.loc[clust_df[f"{rank}"].str.contains("_"), f"{rank}"] = clust_df.loc[
        clust_df[f"{rank}"].str.contains("_"), f"{rank}"
    ] + clust_df.loc[clust_df[f"{rank}"].str.contains("_"), "rank_idx"].astype(str)

## Count sequences in each cluster
clust_df["seq_counts"] = clust_df.groupby("clust_id")["seq_id"].transform("count")

# Split table into single and multiple hits
single_df = clust_df[clust_df.seq_counts == 1]
multi_df = clust_df[clust_df.seq_counts > 1]


# Deal with discrepant clusters (contain multiple taxonomies)
## List unique taxa by cluster for all ranks
multi_df = multi_df.groupby("clust_id", as_index=False)[ranks].aggregate(set)


# Flag clusters with multiple taxa as descrepent for all ranks
## print only 2 species and 4 genus (by default) in formatted_rank
## print all ranks in raw_rank
for rank in ranks:
    multi_df.loc[:, f"formatted_{rank}"] = multi_df.loc[:, f"{rank}"].apply(
        lambda x: format_discrepant_tax(rank, x, ranks_lim)
    )
    multi_df.loc[:, f"all_{rank}"] = multi_df.loc[:, f"{rank}"].apply(
        lambda x: format_discrepant_tax(rank, x, print_desc=False)
    )


## Add seq_id, clust_id, clust_rep and seq_counts to dataframe
multi_df = multi_df.merge(
    clust_df[["clust_id", "seq_id", "clust_rep", "seq_counts"]], on="clust_id"
)

## Columns in raw table
cols = ["clust_id", "clust_rep", "seq_id"] + ranks + ["seq_counts"]


# Df with collapsed taxonomy

formatted_ranks = [f"formatted_{rank}" for rank in ranks]
multi_formatted_df = multi_df[
    ["clust_id", "clust_rep", "seq_id"] + formatted_ranks + ["seq_counts"]
]
multi_formatted_df.columns = cols
collapsed_df = pd.concat([single_df[cols], multi_formatted_df])
collapsed_df["taxpath"] = collapsed_df[ranks].T.agg(";".join)

# Df with uncollapsed taxonomy

all_ranks = [f"all_{rank}" for rank in ranks]
multi_all_df = multi_df[
    ["clust_id", "clust_rep", "seq_id"] + all_ranks + ["seq_counts"]
]
multi_all_df.columns = cols
all_df = pd.concat([single_df[cols], multi_all_df])
all_df["taxpath"] = all_df[ranks].T.agg(";".join)

## Flag discrepant taxa as problematic
all_df.loc[:, "problematic_taxa"] = all_df.apply(problematic_taxa, axis=1)


# Save tables

all_df[all_df.problematic_taxa == True].to_csv(
    snakemake.output.problematic, sep="\t", index=False
)

collapsed_df[["seq_id", "taxpath"]].to_csv(
    snakemake.output.formatted_tax, sep="\t", index=False, header=False
)

all_df[["seq_id", "taxpath"]].to_csv(
    snakemake.output.all, sep="\t", index=False, header=False
)