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

def feat_table(record, p=False):
    """
    Get the feature table from a record as a dictionary.
    If p is set to True, will also print the feature table.
    """
    feat_table = record[0]['GBSeq_feature-table'][0]["GBFeature_quals"]
    d = {}
    for i in feat_table:
        d[i["GBQualifier_name"]] = i["GBQualifier_value"]
    d['tax ID'] = d['db_xref'].split(':')[-1]
    del d['db_xref']
    if p:
        print('\n'.join(f'{k}\t{v}' for k, v in d.items()))
    else:
        return d
