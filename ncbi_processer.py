from Bio import Entrez

"""
This will process the output of our crawler with the Entrez module.
We want the sequence itself, the scientific name, and the tax id.
"""

# ID should be replaced for each
handle = Entrez.efetch(db='nuccore', id='D50540', retmode='xml')
record = Entrez.read(handle)

# This prints a summary of the information
def summary(record):
    print('\n'.join(f'{k}\t{v}' for k, v in record[0].items()))

def feat_table(record):
    feat_table = record[0]['GBSeq_feature-table'][0]["GBFeature_quals"]
    for d in feat_table:
        print('\n'.join(f'{k}\t{v}' for k, v in d.items()))
    return feat_table


table = feat_table(record)

d = {}
for i in table:
    d[i["GBQualifier_name"]] = i["GBQualifier_value"]
