import re, json
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from Bio.SeqIO import parse

#TODO: write parsers for external output from islandviewer4, TAFinder, vfdb

#convience functions for parsing
def copy(handle):
    return handle.read()

def feature_location(position, orientation, contig, start_offset=-1, end_offset=0):
    '''
    Return biopython feature object, where the start position is always lower than the end position, otherwise flip orientation.
    In case of gff files, the strand has to often be overruled.
    :param position: list with two integers [start, end]
    :param no_strand: use True if the direction is not meaningful
    :param strand: overrule the strand, either -1 or 1
    :return: return tuple of (FeatureLocation, strand)
    '''
    f,r = int(position[0]), int(position[1])

    if f <= r:
        location, strand = FeatureLocation(f+start_offset, r+end_offset), 1
    elif f > r:
        location, strand = FeatureLocation(r+start_offset, f+end_offset), -1

    if orientation != None:
        strand = orientation
    return (location, strand, contig)

def primer_length(sequence: str) -> int:
    '''
    Count lenght of primer with square brackets ambigous notation returned by primersearch.
    :param sequence: primer sequence example ATG[TGC]GG?[AT]
    :return: primer lenght
    '''
    # subsitute all square brackets and bases between them for arbitary character
    degenerate = re.sub('\[[A-Z]*\]', '?', sequence)
    return len(degenerate)

def parse_table(handle, header, columns, sep='\t'):
    '''
    Simple table parser. Pandas not used for simplicity and to enable parsing  concatenated files with repeating headers.
    # return list where item represent one line in tab

    :param handle:
    :param header:
    :param columns:
    :param sep:
    :return:
    '''

    return [dict(zip(columns, line.split(sep))) for line in handle if not line.startswith(header)]

def make_records(lines, type_, position, qualifiers, source):
    # create the SeqFeature object
    records = []

    for line in lines:
        location, strand, contig = position(line)
        record = SeqRecord('', contig)
        feature = SeqFeature(FeatureLocation(location[0],location[1], strand=strand), type=type_(line), qualifiers=qualifiers(line))
        record.features = [feature]
        records.append(record)

    return records

#TODO: iterates several times, not particulary effective
def fetch_location(ID, genbank, qualifier='locus_tag'):
    records = parse(genbank, 'genbank')
    location = [(feature.location, feature.strand, record.id)
                for record in records
                for feature in record.features
                if qualifier in feature.qualifiers
                if ID in feature.qualifiers[qualifier]
                ]
    return location[0]

def parse_amplicon(line, primer):
    keys = ('amplicon', 'contig', 'forward_primer','forward_position','forward_mismatches',
            'reverse_primer','reverse_position','reverse_mismatches',  'length'
            )

    x = line.split()
    values = (x[0],x[2],x[3],int(x[8]),int(x[10]),x[12],int(x[17].strip('[]')),int(x[19]),int(x[23]))
    dictionary = dict(zip(keys,values))
    dictionary['primer'] = primer

    return dictionary

#table based outputs are handled by template and helper functions
#TODO: create single function with dictionary of parameters as input
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

def parse_islandpath(handle):
    '''
    Parse islandpath output table file to Biopython Seqfeature objects
    See https://www.pathogenomics.sfu.ca/islandpath/ for details.

    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'islandpath'
    header = '#' #actualy no header here
    colnames = ('name', 'start', 'end')
    type = lambda x: 'genome_island'
    qualifiers = lambda x: {
        'name': x['name'],
        'source': source,
    }
    position = lambda x: feature_location([x['start'], x['end']], 0, x['name'].split('_')[1])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)
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
def parse_REfinder(handle):
    '''
    Parse REfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'REfinder'
    header = 'Gene'
    colnames = (
        'Gene', 'Identity', 'Query_HSP', 'Contig', 'Position_in_contig',
        'Type', 'Function', 'Recognition_Seq','Accession_number'
    )
    type = lambda x: x['Gene']
    qualifiers = lambda x: {
        'name': x['Gene'],
        'score': x['Identity'],
        'alignment': x['Alignment_Length_Gene_Length'],
        'coverage': x['Query_HSP'],
        'reference_location': x['Position_in_reference'],
        'type': x['Type'],
        'product': x['Function'],
        'recognition_seq': x['Recognition_Seq'],
        'db_xref': x['Accession_number'],
        'source': source,
    }
    position = lambda x: feature_location(['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)
def parse_sccmecfinder(handle):
    '''
    Parse sccmecfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'sccmecfinder'
    header = 'Fasta ' #actualy no header here
    colnames = ('Fasta_header', 'Identity', 'Query_HSP_Length', 'Contig', 'Position_in_contig')
    type = lambda x: 'sccmec_cassette'
    qualifiers = lambda x: {
        'name': x['Fasta_header'],
        'score': x['Identity'],
        'alignment': x['Query_HSP_Length'],
        'source': source,
    }
    position = lambda x: feature_location(['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)
def parse_spaTyper(handle):
    '''
    Parse spaTyper output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'spaTyper'
    header = '#' #actualy no header here
    colnames = ('spa_Type', 'Repeats', 'Contig', 'Position', 'Orientation')
    type = lambda x: 'spaType'
    qualifiers = lambda x: {
        'name': x['spa_Type'],
        'repeats': x['Repeats'],
        'source': source,
    }
    position = lambda x: feature_location(['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)
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

#json based outputs are a bit tricky, because of arbitary nesting depth
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
def parse_virulencefinder(handle):
    '''
    Parse virulencefinder json output to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''


    records = []
    results = json.load(handle)['virulencefinder']['results']
    for domain in results:
        for hit in results[domain]:
                if "No hit found" not in results[domain][hit]:
                    for ID in results[domain][hit]:

                        match = results[domain][hit][ID]
                        location, strand, contig = feature_location(match['positions_in_contig'].split('..'),None,match['contig_name'])

                        record = SeqRecord('', contig)
                        qualifiers = {'plasmid': match['plasmid'], #probably need to change this
                                        'identity': match['identity'],
                                        'HSP_length': match['HSP_length'],
                                        'template_length': match['template_length'],
                                        'position_in_ref': match['position_in_ref'],
                                        'note': match['note'],
                                        'db_xref': match['accession'],
                                        'coverage': match['coverage'],
                                        'source':'virulencefinder'
                                      }


                        features = SeqFeature(location, type="gene", strand=strand, qualifiers=qualifiers)
                        record.features = [features]
                        records.append(record)
    return records

#text based inputs are messy and require line by line processing with a lot of exceptions
def parse_primersearch(handle):
    '''
    Parse text file returned by primersearch (Whitehead primer3_core format) to seqrecord object.
    See http://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html for details.
    :param in_file: input file handle
    :return: biopython records list
    '''
    #TODO:separate primer and excluded sequences and create four annotation to each hit

    source = 'primersearch'
    lines = []
    for record in handle.read().strip().split('Primer name '):
        primer = record.split('\nAmplimer ')
        amplicons = [parse_amplicon(amplicon,primer[0])
                      for amplicon in primer[1:]
                      if parse_amplicon(amplicon,primer[0])['length'] <= 3000
                      ]
        lines = lines + amplicons

    records = []
    #amplicons
    type = lambda x: 'amplicon'
    qualifiers = lambda x: {
        'primer': x['primer'],
        'amplicon': x['amplicon'],
        'length': x['length'],
        'mismatches': x['forward_mismatches'] + x['reverse_mismatches'],
        'forward_primer_length': primer_length(x['forward_primer']),
        'reverse_primer_length': primer_length(x['reverse_primer']),
        'forward_primer': x['forward_primer'],
        'reverse_primer': x['reverse_primer'],
        'source': source,
    }
    position = lambda x: feature_location([x['forward_position'], x['forward_position']+x['length']], 0, x['contig'], end_offset=-1)
    records = records + make_records(lines, type, position, qualifiers, source)

    #forward primers
    type = lambda x: 'forward_primer'
    qualifiers = lambda x: {
        'primer': x['primer'],
        'amplicon': x['amplicon'],
        'mismatches': x['forward_mismatches'],
        'length': primer_length(x['forward_primer']),
        'sequence': x['forward_primer'],
        'source': source,
    }
    position = lambda x: feature_location([x['forward_position'], x['forward_position']+primer_length(x['forward_primer'])], 1, x['contig'], end_offset=-1)
    records = records + make_records(lines, type, position, qualifiers, source)

    #reverse primers
    type = lambda x: 'reverse_primer'
    qualifiers = lambda x: {
        'primer': x['primer'],
        'amplicon': x['amplicon'],
        'mismatches': x['reverse_mismatches'],
        'length': primer_length(x['reverse_primer']),
        'sequence': x['reverse_primer'],
        'source': source,
    }
    position = lambda x: feature_location([x['forward_position']+x['length']-primer_length(x['reverse_primer']), x['forward_position']+x['length']], -1, x['contig'], end_offset=-1)
    records = records + make_records(lines, type, position, qualifiers, source)

    #amplicon sequence without primers
    type = lambda x: 'trimmed_amplicon'
    qualifiers = lambda x: {
        'primer': x['primer'],
        'amplicon': x['amplicon'],
        'length': x['length'],
        'mismatches': x['forward_mismatches'] + x['reverse_mismatches'],
        'forward_primer_length': primer_length(x['forward_primer']),
        'reverse_primer_length': primer_length(x['reverse_primer']),
        'forward_primer': x['forward_primer'],
        'reverse_primer': x['reverse_primer'],
        'source': source,
    }
    position = lambda x: feature_location([x['forward_position']+primer_length(x['forward_primer']), x['forward_position']+x['length']-primer_length(x['reverse_primer'])], 0, x['contig'], end_offset=-1)
    records = records + make_records(lines, type, position, qualifiers, source)

    return records

#web based inputs also messy
def parse_web_mge(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'mgefinder'
    header = 'mge_no'
    colnames = (
        'mge_no', 'name', 'synonyms', 'prediction_method', 'type',
        'allelle_len', 'depth', 'e_value', 'identity', 'coverage', 'gaps', 'substitution',
        'contig', 'start', 'stop', 'cigar'
    )
    type = lambda x: x['type']
    qualifiers = lambda x: {
        'name': x['name'],
        'synonyms': x['synonyms'],
        'score': x['e_value'],
        'prediction_method': x['prediction_method'],
        'coverage': x['coverage'],
        'allelle_len': x['allelle_len'],
        'depth': x['depth'],
        'identity': x['identity'],
        'gaps': x['gaps'],
        'substitution': x['substitution'],
        'db_xref': x['cigar'],
        'source': source
    }
    position = lambda x: feature_location([x['start'], x['stop']],None,x['contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_sourcefinder(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'sourcefinder'
    header = '#'
    colnames = (
        'Entry', 'Phage', 'Plasmid', 'Chromosome',
    )
    type = lambda x: x[['Phage', 'Plasmid', 'Chromosome']].idmax(axis=1)
    qualifiers = lambda x: {
        'Phage': x['Phage'],
        'Plasmid': x['Plasmid'],
        'Chromosome': x['Chromosome'],
    }
    position = lambda x: feature_location([[0, -1]],None,x['contig'].strip('>'))
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_plasmidfinder(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'plasmidfinder'
    header = 'Database'
    colnames = (
        'Database', 'Plasmid', 'Identity', 'Query_Template_length',
        'Contig', 'Position_in_contig', 'Note', 'Accession_number',
    )
    type = lambda x: x['Plasmid']
    qualifiers = lambda x: {
        'name': x['Plasmid'],
        'database': x['Database'],
        'identity': x['Identity'],
        'coverage': x['Query_Template'],
        'note': x['Note'],
        'db_xref': x['Accession_number'],
    }
    position = lambda x: feature_location(x['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_sccmecfinder(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'sccmecfinder'
    header = 'Fasta'
    colnames = (
        'Fasta_header', 'Identity', 'Query_HSP_Length', 'Contig','Position_in_contig',
    )
    type = lambda x: x['Plasmid'].split(':')[0]
    qualifiers = lambda x: {
        'name': (':').join(x['Fasta_header'].split(':')[:1]),
        'identity': x['Identity'],
        'coverage': x['Query_HSP_Length'],
        'note': x['Note'],
        'db_xref': (',').join(x['Fasta_header'].split(':')[2:]),
    }
    position = lambda x: feature_location(x['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_refinder(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'sccmecfinder'
    header = 'Fasta'
    colnames = (
        'Gene', 'Identity', 'Query_HSP', 'Contig', 'Position_in_contig'
        'Type', 'Function','Recognition_Seq', 'Accession_number'
    )
    type = lambda x: x['Type']
    qualifiers = lambda x: {
        'name': x['Gene'],
        'identity': x['Identity'],
        'coverage': x['Query_HSP'],
        'product': x['Function'],
        'sequence': x['Recognition_Seq'],
        'db_xref': x['Accession_number'].strip(','),
    }
    position = lambda x: feature_location(x['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_vfdb(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'plasmidfinder'
    header = 'Database'
    colnames = (
        'Database', 'Plasmid', 'Identity', 'Query_Template_length',
        'Contig', 'Position_in_contig', 'Note', 'Accession_number',
    )
    type = lambda x: x['Plasmid']
    qualifiers = lambda x: {
        'name': x['Plasmid'],
        'database': x['Database'],
        'identity': x['Identity'],
        'coverage': x['Query_Template'],
        'note': x['Note'],
        'db_xref': x['Accession_number'],
    }
    position = lambda x: feature_location(x['Position_in_contig'].split('..'),None,x['Contig'])
    lines = parse_table(handle, header, colnames)

    return make_records(lines, type, position, qualifiers, source)

def parse_web_islandviewer(handle):
    '''
    Parse resfinder output table from the web application to Biopython Seqfeature objects
    See  for details.
    :param table: file handle
    :return: list of SeqRecords
    '''

    source = 'islandviewer'
    header = 'Start'
    colnames = (
        'Start', 'End', 'Size', 'GI_Prediction_Method', 'View',
        'External_Annotations', 'Download',
    )
    type = lambda x: x['GI_Prediction_Method']
    qualifiers = lambda x: {
        'db_xref': x['External_Annotations'].strip(),
    }
    position = lambda x: feature_location([x['Start'].replace(',',''),x['End'].replace(',','')],None,None)
    lines = parse_table(handle, header, colnames)
    return make_records(lines, type, position, qualifiers, source)

nucleotide = {
    'abricate': parse_abricate,
    'islandpath':parse_islandpath,
    'rgi':parse_rgi,
    'resfinder':parse_resfinder,
    'refinder':parse_REfinder,
    'sccmecfinder':parse_sccmecfinder,
    'spatyper':parse_spaTyper,
    'plasmidfinder':parse_plasmidfinder,
    'crisprcasfinder':parse_crisprcasfinder,
    'virulencefinder':parse_virulencefinder,
    'primersearch': parse_primersearch,
    'islandviewer': parse_web_islandviewer,
    'mge': parse_web_mge,
    #'source':parse_web_sourcefinder,
}
protein = {
    'bacmet2': parse_bacmetscan,
    'web_vfdb': parse_web_vfdb,
}
copy = {
    #'mge':copy,
    'keyword':copy,
    'phispy':copy,
}

if __name__ == "__main__":
    genbank = snakemake.input[0]
    annotations = snakemake.input[1:]
    for annotation in annotations:
        print('Make', annotation)
        program = annotation.split('.')[-2].split('_')[0]
        with open(annotation, 'r') as in_file:
            print(program)
            print([output for output in snakemake.output])
            print([output for output in snakemake.output if program in output])
            [output for output in snakemake.output if program in output.split('.')[-2]]
            output = [output for output in snakemake.output if program in output.split('.')[-2]][0]
            print(output)
            with open(output, 'w') as out_file:
                if program in copy:
                    out_file.write(copy[program](in_file))
                else:
                    if program in nucleotide:
                        records = nucleotide[program](in_file)
                    elif program in protein:
                        records = protein[program](in_file, genbank)
                    GFF.write(records, out_file)
