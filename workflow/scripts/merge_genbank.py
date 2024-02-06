from BCBio import GFF
from Bio import SeqIO

#TODO: genbank add logic for merging genbank files

def merge_genbank(record, annotation, import_id=True, to_import=['gene'], no_import=['source']):
    '''
    Merge annotations from two genbank files. Files have to contain the same sequence.
    It is intended to merge rast - record and prokka - annotation
    :param record: base record
    :param annotation: record with annotations
    :param import_id: import id boolean, default True
    :param to_import: features to import
    :param no_import: features to not import
    :return: SeqRecord with merged features
    '''

    assert record.seq == annotation.seq #check if both records contain the same sequence

    if import_id:
        record.id = annotation.id

    features_generator = (feature
                          for feature in Annotation.features
                          if feature.type not in no_import
                          if feature.type in to_import
                          )
    for feature in features_generator:
        record.features.append(feature)

    return record

def merge_gff(genbank_file, gff_files):
    '''
    Merge multiple gff files into one genbank file. Files have to contain the same record/contig ids.
    If none provided, its assumed that it belongs to first record.
    It is intended to import multiple gff annotations into prokka gbk file.
    Beware that the sequence identity is not checked (can be problem if flipped)
    :param record: path to genbank file
    :param annotation: list of paths to gff files
    :return: list of SeqRecords with merged features
    '''

    #load all gff records

    gff_records = []
    for gff_file in gff_files:
        handle = open(gff_file, 'r')
        gff = [record for record in GFF.parse(handle)]
        gff_records = gff_records + gff
        handle.close()

    #append features from gff records to genbank records
    new_records = []
    first_record = True
    for record in SeqIO.parse(genbank_file, 'genbank'):
        for gff_record in gff_records:
            if record.id == gff_record.id:
                for feature in gff_record.features:
                    record.features.append(feature)
            elif gff_record.id == 'None' and first_record == True:
                print('triggered')
                for feature in gff_record.features:
                    record.features.append(feature)
            else:
                print(record.id)
        first_record = False

        new_records.append(record)

    return new_records

if __name__ == "__main__":
    SeqIO.write(merge_gff(snakemake.input[0], snakemake.input[1:]), snakemake.output[0], 'genbank')
