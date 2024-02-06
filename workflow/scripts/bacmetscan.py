from BCBio import GFF
from AnyToGff import feature_location, parse_table, make_records, fetch_location


def parse_bacmetscan(handle, genbank):
    '''
    Parse bacmetscan output table file to Biopython Seqfeature objects
    See http://bacmet.biomedicine.gu.se/download_temporary.html for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'bacmetscan'
    header = 'Query'
    colnames = (
        'Query', 'Subject', 'Gene', 'Description',
        'Organism', 'Location', 'Compounds', 'Percent_identity',
        'Match_length', 'E_value', 'Score_per_length'
    )
    type = lambda x: 'protein'
    qualifiers = lambda x: {
        'name': x['Subject'],
        'coverage': x['Match_length'],
        'identity': x['Percent_identity'],
        'db_xref': x['Gene'],
        'resistance_phenotype': x['Compounds'],
        'organism': x['Location'],
        'type': x['Organism'],
        'description': x['Description'],
        'e_value': x['E_value'],
        'score': x['Score_per_length'],
        'source': source,
    }
    position = lambda x: fetch_location(x['Query'], genbank)

    #the bacmetscan has funny shaped output table, possibly bug
    lines = parse_table(handle.read().strip('\n').replace('\n\t', '\t').split('\n'), header, colnames)

    return make_records(lines, type, position, qualifiers, source)

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_bacmetscan(in_file, genbank)
            GFF.write(records, out_file)