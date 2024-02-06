from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from AnyToGff import feature_location, parse_table, make_records
import json

def parse_plasmidfinder(handle):
    '''
    Parse plasmidfinder json output to Biopython Seqfeature objects
    See for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    records = []
    results = json.load(handle)['plasmidfinder']['results']
    for domain in results:
        for hit in results[domain]:
                if "No hit found" not in results[domain][hit]:
                    for ID in results[domain][hit]:

                        match = results[domain][hit][ID]
                        location, strand, contig = feature_location(match['positions_in_contig'].split('..'), None, match['contig_name'])

                        record = SeqRecord('', contig)
                        qualifiers = {'plasmid': match['plasmid'],
                                        'identity': match['identity'],
                                        'HSP_length': match['HSP_length'],
                                        'template_length': match['template_length'],
                                        'position_in_ref': match['position_in_ref'],
                                        'note': match['note'],
                                        'db_xref': match['accession'],
                                        'coverage': match['coverage'],
                                        'source':'plasmidfinder'
                                      }

                        features = SeqFeature(location, type="gene", strand=strand, qualifiers=qualifiers)
                        record.features = [features]
                        records.append(record)
    return records

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_plasmidfinder(in_file)
            GFF.write(records, out_file)