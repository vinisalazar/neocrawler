from Bio import Entrez

"""
This will process the output of our crawler with the Entrez module.
We want the sequence itself, the scientific name, and the tax id.
"""

# ID should be replaced for each
handle = Entrez.efetch(db='nuccore', id='D50540', retmode='xml')
record = Entrez.read(handle)

print('\n'.join(f'{k}\t{v}' for k, v in record[0].items()))
print('\n'.join(f'{k}\t{v}' for k, v in record[0]['GBSeq_feature-table'][0].items()))
