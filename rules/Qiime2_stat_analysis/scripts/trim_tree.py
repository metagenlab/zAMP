import qiime2
from skbio import TreeNode
import pandas as pd
from gneiss.util import match_tips


table = qiime2.Artifact.load(snakemake.input["table"]).view(pd.DataFrame)
tree = qiime2.Artifact.load(snakemake.input["tree"]).view(pd.DataFrame)
table, tree = match_tips(table, tree)
tree.write(snakemake.output)
