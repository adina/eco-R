import sys 
from Bio import SeqIO

#to run on multiple, run
#for x in *gbk; do python parse-genbank2.py $x > $x.16s.fa; done

genome=SeqIO.read(sys.argv[1], 'genbank')

n = 0
l = []
for feat in genome.features:
    if feat.type == "rRNA":
        if feat.qualifiers['product'][0] == '16S ribosomal RNA':
            start = feat.location.start.position
            end = feat.location.end.position
            pos = [start, end]
            l.append(pos)
            print '>' + sys.argv[1].split('.')[0] + ' 16S rRNA gene ' + str(n) 
            print feat.extract(genome.seq)
            n =+ 1

