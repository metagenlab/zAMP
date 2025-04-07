import pandas as pd
from snakemake.script import snakemake


# Functions
def shorten_taxa(tax, n=4):
    if "/" in tax:
        tax_list = tax.split("/")
        if len(tax_list) > n:
            return "/".join(tax_list[0:n]) + f"/({len(tax_list) - n})"
        else:
            return tax
    else:
        return tax


def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)


def format_species(species):
    sp = species.replace("_s", " s")
    unique_genera = list(
        set([name.split(" ")[0].split("_s")[0] for name in sp.split("/")])
    )
    if len(unique_genera) == 1:
        genus = unique_genera[0]
        return rreplace(sp, genus, "", sp.count(genus) - 1).replace("/ ", "/")
    else:
        return species


def ambiguous_taxa(row):
    return any("/" in str(cell) for cell in row)


# Ranks list

ranks = snakemake.params[0].split(",")

# Read tables
## Taxonomy table
tax_df = pd.read_csv(snakemake.input[0], sep="\t")
if any(tax_df.columns.str.contains(";")):
    tax_df = pd.read_csv(snakemake.input[0], sep="\t", header=None)
tax_df.columns = ["seq_id", "tax"]
## Split taxonomy path into one column per rank
tax_df[ranks] = tax_df.tax.str.split(";", expand=True).loc[:, 0:6]


# Vsearch output table
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
clust_df = clust_df.replace("_placeholder", "", regex=True)

## Count sequences in each cluster
clust_df["seq_counts"] = clust_df.groupby("clust_id")["seq_id"].transform("count")

# Split table into single and multiple sequence clusters
single_df = clust_df[clust_df.seq_counts == 1]
# Count taxa occurences in multi-sequence clusters
multi_df = (
    clust_df[clust_df.seq_counts > 1]
    .groupby("clust_id", as_index=False)[ranks]
    .value_counts()
)

# Sort taxa occurence for each cluster
multi_df = multi_df.sort_values(
    ["clust_id", "count", "species"], ascending=[True, False, True]
)

# All taxa names are joined in "/" seperated name
# Most frequent taxa in the cluster appear first
multi_df = (
    multi_df.groupby("clust_id")[ranks]
    .agg(lambda x: "/".join(dict.fromkeys(x)))
    .reset_index()
)

multi_df["formatted_Species"] = multi_df.species.apply(lambda x: format_species(x))
multi_df["short_format_Species"] = multi_df.formatted_Species.apply(
    lambda x: shorten_taxa(x)
)
## Add seq_id to multi_df
multi_df = clust_df[clust_df.record_type == "C"][["clust_id", "seq_id"]].merge(
    multi_df, on="clust_id"
)

all_tax_df = pd.concat([single_df[["seq_id"] + ranks], multi_df[["seq_id"] + ranks]])
all_tax_df.loc[:, "ambiguous_taxa"] = all_tax_df.apply(ambiguous_taxa, axis=1)
all_tax_df["all_tax"] = all_tax_df[ranks].T.agg(";".join)

multi_form_tax_df = multi_df[["seq_id"] + ranks[:-1] + ["short_format_Species"]]
multi_form_tax_df = multi_form_tax_df.rename(
    columns={"short_format_Species": "species"}
)
formatted_tax_df = pd.concat([single_df[["seq_id"] + ranks], multi_form_tax_df])
formatted_tax_df["formatted_tax"] = formatted_tax_df[ranks].T.agg(";".join)

## Sort sequence IDs as they appear in the fasta file
seq_id_order = clust_df[clust_df.record_type == "C"]["seq_id"]
all_tax_df["seq_id"] = pd.Categorical(
    all_tax_df["seq_id"], categories=seq_id_order, ordered=True
)
all_tax_df = all_tax_df.sort_values("seq_id")

formatted_tax_df["seq_id"] = pd.Categorical(
    formatted_tax_df["seq_id"], categories=seq_id_order, ordered=True
)
formatted_tax_df = formatted_tax_df.sort_values("seq_id")


all_tax_df[all_tax_df.ambiguous_taxa].to_csv(
    snakemake.output.ambiguous, sep="\t", index=False
)

formatted_tax_df[["seq_id", "formatted_tax"]].merge(
    all_tax_df[["seq_id", "all_tax"]], on="seq_id"
).to_csv(snakemake.output.all, sep="\t", index=False)


formatted_tax_df[["seq_id", "formatted_tax"]].to_csv(
    snakemake.output.collapsed, sep="\t", index=False, header=False
)
