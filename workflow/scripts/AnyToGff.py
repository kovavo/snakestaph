import re
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
        feature = SeqFeature(FeatureLocation(int(location.start),int(location.end), strand=strand), type=type_(line), qualifiers=qualifiers(line))
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