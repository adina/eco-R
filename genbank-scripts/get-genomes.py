from Bio import Entrez
import sys

l = []
for line in open(sys.argv[1]):
    l.append(line.rstrip())

for each in l:
    Entrez.email = "howead@msu.edu"
    handle = Entrez.efetch(db="nucleotide", id=each, rettype="gb", retmode="text")
    fp = open(each + '.gbk', 'w')
    fp.write('%s' % handle.read())
