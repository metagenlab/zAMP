#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 3:
    exit("Usage: get_f1_scores.py <tsv> <prefix>")


def get_match(expected, assigned):
    return all(word in str(assigned) for word in str(expected).split())


def classify_match(row, rank):
    return "tp" if get_match(row[f"{rank}_expected"], row[f"{rank}_assigned"]) else "fp"


def compute_metrics(group, rank):
    tp = (group[f"{rank}_prediction"] == "tp").sum()
    fp = (group[f"{rank}_prediction"] == "fp").sum()
    fn = int(group[f"{rank}_fn"].sum())
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_score = (
        2 * (precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0
    )

    return pd.Series(
        {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1_score": f1_score,
        }
    )


df = pd.read_csv(sys.argv[1], sep="\t")
df = df[df.amplified]

ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
all_scores = []
for rank in ranks:
    df[f"{rank}_prediction"] = df.apply(
        lambda row: classify_match(row, rank=rank),
        axis=1,
    )
    assigned = df[f"{rank}_assigned"].tolist()
    df[f"{rank}_fn"] = df[f"{rank}_expected"].apply(
        lambda x: not any(get_match(x, assigned) for assigned in assigned)
    )
    scores = df.groupby("fasta", as_index=False).apply(
        compute_metrics, rank, include_groups=False
    )
    scores["rank"] = [rank] * len(scores)
    all_scores.append(scores)
all_scores = pd.concat(all_scores)


all_scores.to_csv(f"{sys.argv[2]}_all_scores.tsv", sep="\t", index=False)

mean_scores = all_scores.groupby("rank", as_index=False).agg(
    {
        "tp": "sum",
        "fp": "sum",
        "fn": "sum",
        "precision": "mean",
        "recall": "mean",
        "f1_score": "mean",
    }
)
mean_scores["rank"] = pd.Categorical(
    mean_scores["rank"],
    categories=ranks,
    ordered=True,
)
mean_scores = mean_scores.sort_values("rank")
mean_scores.to_csv(f"{sys.argv[2]}_mean_scores.tsv", sep="\t", index=False)
