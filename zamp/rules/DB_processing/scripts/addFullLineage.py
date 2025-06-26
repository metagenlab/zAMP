## Reproduced from https://github.com/GLBRC-TeamMicrobiome/python_scripts

#!/usr/bin/python
import sys

if len(sys.argv) != 3:
    exit("addFullLineage.py taxonomyFile fastaFile")

f1 = open(sys.argv[1], "r").readlines()
hash = {}  # lineage map
for line in f1[1:]:
    line = line.strip()
    # convert to unicode to avoid trouble from non-unicode source file
    # line = unicode(line, "UTF-8")#strip non-unicode non-breaking space
    # line = line.strip("\u00a0")
    cols = line.strip().split("\t")
    lineage = ["Root"]
    for node in cols[1:]:
        node = node.strip()
        if not (node == "-" or node == ""):
            lineage.append(node)
    ID = cols[0]
    lineage = ";".join(lineage).strip()
    hash[ID] = lineage
f2 = open(sys.argv[2], "r").readlines()
for line in f2:
    line = line.strip()
    if line == "":
        continue
    # convert to unicode to avoid trouble from non-unicode source file
    # line = unicode(line, "UTF-8")#strip non-unicode non-breaking space
    # line = line.strip(u"\u00A0")
    if line[0] == ">":
        ID = line.strip().split()[0].replace(">", "")
        try:
            lineage = hash[ID]
        except KeyError:
            print(ID, "not in taxonomy file")
        print(">" + ID + "\t" + lineage)
    else:
        print(line.strip())
