from Bio.SeqIO import parse
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

def keyword2genbank(genbank, features,  qualifiers, values):

    #parsing genbank file, maybe should be a part of different function
    records = parse(genbank, "genbank")

    for record in records:
        #list each feature which contain keywords in given qualifiers

        key_features = [(feature, keyword, value)
                        for feature in record.features
                        if feature.type in features
                        for qualifier in qualifiers
                        for keyword in values
                        if qualifier in feature.qualifiers
                        for value in feature.qualifiers[qualifier]
                        if keyword.lower() in value.lower()
                        ]

        new_record = SeqRecord('', record.id)

        #create new annotation and add it to the record
        for feature, keyword, value in key_features:
            new_feature = SeqFeature(feature.location, type=keyword, qualifiers={'name': keyword, 'source':'keyword'})
            new_record.features.append(new_feature)

        yield new_record

with open(snakemake.output[0], "w") as out_handle:
    GFF.write(keyword2genbank(snakemake.input[0], snakemake.params[0], snakemake.params[1], snakemake.params[2]), out_handle)
