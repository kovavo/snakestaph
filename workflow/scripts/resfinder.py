from BCBio import GFF
from AnyToGff import feature_location, parse_table, make_records

def parse_resfinder(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'resfinder'
    header = 'Resistance'
    colnames = (
        'Resistance_gene', 'Identity', 'Alignment_Length_Gene_Length', 'Coverage',
        'Position_in_reference', 'Contig', 'Position_in_contig', 'Phenotype', 'Accession'
    )
    type = lambda x: x['Resistance_gene']
    qualifiers = lambda x: {
        'name': x['Resistance_gene'],
        'score': x['Identity'],
        'alignment': x['Alignment_Length_Gene_Length'],
        'coverage': x['Coverage'],
        'reference_location': x['Position_in_reference'],
        'resistance_phenotype': x['Phenotype'],
        'db_xref': x['Accession'],
        'source': source,
    }
    position = lambda x: feature_location(x['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_resfinder(in_file)
            GFF.write(records, out_file)