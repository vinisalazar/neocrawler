from Bio import SeqIO

"""
This script dereps a .fasta file.
TODO: add command line functionality for input and output parameters.
"""

records = list(SeqIO.parse("data/seqs_formatted_description.fasta", "fasta"))
enum = list(enumerate(records))

d = {}
for n in enum:
     d[n[0]] = n[1]

drop_d = {}

for key,value in d.items():
    if value.description not in [i.description for i in drop_d.values()]:
        drop_d[key] = value

records = []
for v in drop_d.values():
    records.append(v)


with open("seqs_derep_1.fasta", "w") as f:
    SeqIO.write(records, f, "fasta")
