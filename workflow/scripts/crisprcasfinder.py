from BCBio import GFF
from AnyToGff import feature_location, parse_table, make_records
import json


def parse_crisprcasfinder(handle):
    '''
    Parse crisprcasfinder++ json output to SeqFeature objects
    See https://crisprcas.i2bc.paris-saclay.fr/ for details.
    :param handle: file handle
    :return: list of SeqRecords
    '''

    records = []
    #orientation is either +,- or ND, if ND return 0
    orient = lambda x: x if x not in 'ND' else '0'

    for sequence in json.load(handle)['Sequences']:
        #create empty record object for each sequence
        source = 'crisprcasfinder'
        contig = sequence['Id']


        #Cas
        lines = sequence['Cas']
        type = lambda x: 'cas_region'
        qualifiers = lambda x: {
            'type':x['Type'],
            'source': source
        }
        position = lambda x: feature_location([x['Start'], x['End']], 0, contig)

        cas_regions = make_records(lines, type, position, qualifiers, source)

        #Cas genes
        lines = [gene for cas in sequence['Cas'] for gene in cas['Genes']]
        type = lambda x: 'Cas'
        qualifiers = lambda x: {
            'type':x['Sub_type'],
            'source': source
        }
        position = lambda x: feature_location([x['Start'], x['End']], int(orient(x['Orientation'])+'1'), contig)

        cas_genes = make_records(lines, type, position, qualifiers, source)

        #CRIPSR region
        lines = sequence['Crisprs']
        type = lambda x: 'crispr_region'
        qualifiers = lambda x: {
            'id':x['Name'],
            'DR_conservation':x['Conservation_DRs'],
            'DR_consensus':x['DR_Consensus'],
            'DR_length': x['DR_Length'],
            'evidence_level': x['Evidence_Level'],
            'spacers': x['Spacers'],
            'repeat_ID': x['Repeat_ID'],
            'source': source
        }
        position = lambda x: feature_location([x['Start'], x['End']], int(orient(x['Potential_Orientation'])+'1'), contig)

        crispr_regions = make_records(lines, type, position, qualifiers, source)

        #CRISPR elements
        lines = [region for crispr in sequence['Crisprs'] for region in crispr['Regions']]
        type = lambda x: x['Type']
        qualifiers = lambda x: {
            'sequence':x['Sequence'],
            'source': source
        }
        position = lambda x: feature_location([x['Start'], x['End']], 0, contig)

        crispr_elements = make_records(lines, type, position, qualifiers, source)

        records = records + cas_regions + cas_genes + crispr_regions + crispr_elements

    return records

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_crisprcasfinder(in_file)
            GFF.write(records, out_file)