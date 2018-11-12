from Bio import SeqIO

records = list(SeqIO.parse("seqs_derep_trimmed_dedup.fasta", "fasta"))

ix = 000000
for seq in records:
    seq.id = str(ix)
    ix += 1
    seq.name = f"{seq.id} " + seq.description.split(',')[0]
    seq.dbxrefs = int(seq.description.split()[-1])

with open("bacterio.map", "w") as f:
    for seq in records:
        f.write(f"{seq.name}\t{seq.dbxrefs}\n")

with open("bacterio.fasta", "w") as f:
    SeqIO.write(records, f, "fasta")
