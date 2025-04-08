#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 4:
    exit("Usage: format_RDP_output.py <ranks> <input> <output>")

ranks = sys.argv[1].split(",")
df = pd.read_csv(sys.argv[2], sep="\t", header=None)
ranks_idx = list(range(5, (len(ranks) - 1) * 4, 3))
conf_idx = ranks_idx[-1] + 2
cols = [0] + ranks_idx + [conf_idx]
df = df[cols]
df.columns = ["seq_id"] + ranks + ["confidence"]
df.insert(1, "tax", df[ranks].T.agg(";".join))
df[["seq_id", "tax", "confidence"]].to_csv(
    sys.argv[3], sep="\t", index=False, header=False
)
