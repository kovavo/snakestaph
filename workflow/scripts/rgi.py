from BCBio import GFF
from AnyToGff import feature_location, parse_table, make_records

def parse_rgi(handle):
    '''
    Parse rgi/card output json file to Biopython Seqfeature objects
    See https://card.mcmaster.ca/analyze/rgi for details.

    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'rgi'
    header = 'ORF'
    colnames = (
        'ORF_ID', 'Contig', 'Start', 'Stop', 'Orientation', 'Cut_Off', 'Pass_Bitscore', 'Best_Hit_Bitscore',
        'Best_Hit_ARO', 'Best_Identities', 'ARO', 'Model_type', 'SNPs_in_Best_Hit_ARO', 'Other_SNPs', 'Drug_Class',
        'Resistance_Mechanism', 'AMR_Gene_Family', 'Predicted_DNA', 'Predicted_Protein', 'CARD_Protein_Sequence',
        'Percentage_Length_of_Reference_Sequence', 'ID', 'Model_ID', 'Nudged', 'Note'
    )
    type = lambda x: x['Cut_Off']
    qualifiers = lambda x: {
        'name': x['Best_Hit_ARO'],
        'identity': x['Best_Identities'],
        'db_xref_ARO': x['ARO'],
        'pass_bitscore': x['Pass_Bitscore'],
        'best_bitscore': x['Best_Hit_Bitscore'],
        'model_type': x['Model_type'],
        'best_snps': x['SNPs_in_Best_Hit_ARO'],
        'other_snps': x['Other_SNPs'],
        'drug_class': x['Drug_Class'],
        'drug_mechanism': x['Resistance_Mechanism'],
        'product': x['AMR_Gene_Family'],
        'predicted_DNA': x['Predicted_DNA'],
        'predicted_protein': x['Predicted_Protein'],
        'CARD_protein': x['CARD_Protein_Sequence'],
        'coverage': x['Percentage_Length_of_Reference_Sequence'],
        'db_xref': x['ID'],
        'model_ID': x['Model_ID'],
        'nudged': x['Nudged'],
        'note': x['Note'],
        'source': source,
    }
    position = lambda x: feature_location([x['Start'], x['Stop']], int(x['Orientation']+'1'),'_'.join(x['Contig'].split('_')[0:-1]))
    lines = parse_table(handle, header, colnames)
    return make_records(lines, type, position, qualifiers, source)

if __name__ == "__main__":
    genbank, annotation, output = snakemake.input[0], snakemake.input[1], snakemake.output[0]
    print('Converting', annotation, 'output to gff')
    with open(annotation, 'r') as in_file:
        with open(output, 'w') as out_file:
            records = parse_rgi(in_file)
            GFF.write(records, out_file)