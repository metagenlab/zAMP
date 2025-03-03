#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 3:
    exit("Usage: format_RDP_output.py <input> <output>")

df = pd.read_csv(sys.argv[1], sep="\t", header=None)
ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

out_df = df[[0, 5, 8, 11, 14, 17, 20, 23, 25]]
out_df.columns = ["seq_id"] + ranks + ["confidence"]
out_df.insert(1, "tax", out_df[ranks].T.agg(";".join))
out_df[["seq_id", "tax", "confidence"]].to_csv(
    sys.argv[2], sep="\t", index=False, header=False
)
