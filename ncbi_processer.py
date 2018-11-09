__author__ = "Vini Salazar"
import os
import argparse
from Bio import Entrez

# Please replace with your own e-mail
Entrez.email = "vinicius.salazar@neoprospecta.com"

"""
This will process the output of our crawler with the Entrez module.
    - cleaner(file) parses the crawler output for downstream processing.
    - summary(record) prints a summary of a record contained in our file.
    - feat_table(record) returns the feature table for that record.
    - pipeline(kind) executes the output processing pipeline. 'kind' specifies
      how we want that to be done.

"""

def cleaner(file, overwrite=True, out=None):
    """
    This function cleans our crawler output. If overwritten is false, the
    original output will be kept. Otherwise, it will be cleaned. Use the out
    kwarg to specify a suffix to the outputfile.
    """
    with open(file, 'r') as f:
        r = f.readlines()[1:]
        r = [i.strip() for i in r if len(i) > 4 and not len(i.split()) > 1]
        r = [i.replace('"', '') for i in r]
        r = [i.replace(',', '') for i in r]

    if out:
        file = file.split('.')[-2] + out + '.' + file.split('.')[-1]

    if overwrite:
        with open(file, 'w') as f:
            [f.write(f'{i}\n') for i in r]

    return r


def ncbi_request(record, p=False, feats=False):
    """
    Requests a record (accession number) from the NCBI Entrez and reads it.
    If summary is passed as True, prints the record with summary()
    If feat_table is passed as True, prints the feature table with feat_table()
    """
    handle = Entrez.efetch(db='nuccore', id=record, format='xml')
    record = Entrez.read(handle)

    if p:
        summary(record)

    if feats:
        feat_table(record, p=True)

    return record


# This prints a summary of the information
def summary(record):
    """
    Prints a summary of the record associated with our accession number.
    'references', 'feature_table' and 'sequence' are excluded because they
    are too 'noisy'.
    """
    print('\n'.join(f'{k}\t{v}'
          for k, v in record[0].items() if
              'GBSeq_references' not in k and
              'GBSeq_feature-table' not in k and
              'GBSeq_sequence' not in k)
    )


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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script processes our crawler output."
        )
    parser.add_argument("-i", "--input", type=str, default="neocrawler.csv")
    args = parser.parse_args()

    records = cleaner(args.input, overwrite=False)
    for record in records:
        ncbi_request(record, p=True, feats=True)
