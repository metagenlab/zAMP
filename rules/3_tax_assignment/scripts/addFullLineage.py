#!/usr/bin/python
import sys, string
if len(sys.argv) != 3:
	print 'addFullLineage.py taxonomyFile fastaFile'
	sys.exit()



input_file_1 = snakemake.input[0]
input_file_2 = snakemake.input[1]
output_file = snakemake.output[0]



f1 = open(input_file_1, 'r').readlines()
hash = {} #lineage map
for line in f1[1:]:
	line = line.strip()
	#convert to unicode to avoid trouble from non-unicode source file
	#line = unicode(line, "UTF-8")#strip non-unicode non-breaking space
	#line = line.strip(u"\u00A0")
	cols = line.strip().split('\t')
	lineage = ['Root']
	for node in cols[1:]:
		node = node.strip()
		if not (node == '-' or node == ''):
			lineage.append(node)
	ID = cols[0]
	lineage = string.join(lineage, ';').strip()
	hash[ID] = lineage
f2 = open(input_file_2, 'r').readlines()

f = open(output_file, 'w')

saveout = sys.stdout 
sys.stdout = f

for line in f2:
	line = line.strip()
	if line == '':
		continue
	#convert to unicode to avoid trouble from non-unicode source file
	#line = unicode(line, "UTF-8")#strip non-unicode non-breaking space
	#line = line.strip(u"\u00A0")
	if line[0] == '>':
		ID = line.strip().split()[0].replace('>', '')
		try:
			lineage = hash[ID]
		except KeyError:
			print ID, 'not in taxonomy file'
			sys.exit()
		print '>' + ID + '\t' + lineage
	else:
			print line.strip()

sys.stdout = saveout 
f.close()      