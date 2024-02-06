from BCBio import GFF
from AnyToGff import feature_location, parse_table, make_records

def parse_abricate(handle):
    '''
    Parse abricate output table file to Biopython Seqfeature objects
    see https://github.com/tseemann/abricate for details

    :param table: file handle
    :return: list of SeqRecords
    '''
    source='abricate'
    header = '#'
    colnames = (
        'file','sequence', 'start', 'end', 'strand',
        'gene', 'count', 'map', 'gaps', 'coverage',
        'identity', 'database', 'accession', 'product', 'resistance'
    )
    lines = parse_table(handle, header, colnames)
    type = lambda x: x['database']
    qualifiers = lambda x:{
        'name': x['gene'],
        'coverage': x['coverage'],
        'identity': x['identity'],
        'db_xref': x['accession'],
        'product': x['product'],
        'resistance_phenotype':x['resistance'],
        'coverage_map':x['map'],
        'coverage_count':x['count'],
        'gaps_count':x['gaps'],
        'source':source,
    }
    position = lambda x: feature_location([x['start'], x['end']], int(x['strand']+'1'), x['sequence'])

    return make_records(lines, type, position, qualifiers, source)

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_abricate(in_file)
            GFF.write(records, out_file)