__author__ = "Vini Salazar"
import os
import argparse
from time import sleep
from Bio import SeqIO
from Bio import Entrez

# Please replace with your own e-mail
Entrez.email = "myemail@email.com"

"""
This will process the output of our crawler with the Entrez module.
    - cleaner(file) parses the crawler output for downstream processing.
    - summary(record) prints a summary of a record contained in our file.
    - feat_table(record) returns the feature table for that record.
    - pipeline(kind) executes the output processing pipeline. 'kind' specifies
      how we want that to be done.

"""

def cleaner(file, overwrite=False, out=None, skip=True):
    """
    This function cleans our crawler output. If overwritten is false, the
    original output will be kept. Otherwise, it will be cleaned. Use the out
    kwarg to specify a suffix to the outputfile.
    """
    with open(file, 'r') as f:
        if skip:
            r = f.readlines()[1:]
        else:
            r = f.readlines()

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

    return record[0]


# This prints a summary of the information
def summary(record):
    """
    Prints a summary of the record associated with our accession number.
    'references', 'feature_table' and 'sequence' are excluded because they
    are too 'noisy'.
    """
    summary = ('\n'.join(f'{k}\t{v}'
          for k, v in record.items() if
              'GBSeq_references' not in k and
              'GBSeq_feature-table' not in k and
              'GBSeq_sequence' not in k)
    )

    return summary


def feat_table(record, p=False):
    """
    Get the feature table from a record as a dictionary.
    If p is set to True, will also print the feature table.
    """
    feat_table = record['GBSeq_feature-table'][0]["GBFeature_quals"]
    d = {}
    for i in feat_table:
        d[i["GBQualifier_name"]] = i["GBQualifier_value"]
    d['tax ID'] = d['db_xref'].split(':')[-1]
    del d['db_xref']
    if p:
        print('\n'.join(f'{k}\t{v}' for k, v in d.items()))
    else:
        return d


def writer(record, taxfile, sequencefile, recordsfile):

    taxfile, sequencefile, recordsfile = os.path.join('data/', taxfile),\
                                         os.path.join('data/', sequencefile),\
                                         os.path.join('data/', recordsfile)

    features = feat_table(record)
    info = summary(record)
    try:
        name, strain, taxid = features['organism'],\
                            features['strain'],\
                            features['tax ID']
    except KeyError:
        name, strain, taxid = features['organism'],\
                              '<Unknown strain>',\
                              features['tax ID']

    header = str(features)[13:-2].replace("'", '')
    try:
        seq = record["GBSeq_sequence"].upper()
    except KeyError:
        seq = 'ATCG'

    with open(taxfile, 'a') as f:
        f.write(f'{name} strain {strain}\t{taxid}\n')
        f.close()

    with open(sequencefile, 'a') as f:
        f.write(f'> {header}\n{seq}\n')
        f.close()

    with open(recordsfile, 'a') as f:
        f.write(info + '\n')
        f.close()

    return f'Written {record} to {taxfile}, {sequencefile}, {recordsfile}'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script processes our crawler output."
        )

    parser.add_argument("-i", "--input", type=str, default="neocrawler.csv",
                        help="Path to input file in .csv format.\
                             can be either clean or dirty crawler output.")
    parser.add_argument("-r", "--recordsfile", default="raw_records.csv",
                        help="Path to output file with raw records.")
    parser.add_argument("-t", "--taxfile", type=str, default="crawler_tax.csv",
                        help="Path to output tax ID and tax names file.")
    parser.add_argument("-s", "--sequencefile", type=str, default="crawler_seqs.fasta",
                        help="Path to output fasta file containing sequences.")

    args = parser.parse_args()

    records = cleaner(args.input)

    for record in records:
        try:
            record = ncbi_request(record)
            print(summary(record))
            writer(record, args.taxfile, args.sequencefile, args.recordsfile)
            sleep(2)
        except:
            pass
