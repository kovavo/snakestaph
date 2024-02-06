'''
Master snakefile for all workflows
'''

import pandas as pd
import os

config['paralel'] = 1

def targets(df, mode):

    '''
    Create list of all target files based on samples and chosen pipeline
    '''

    ##################################
    ##        reads as input        ##
    ##################################

    if mode in ['quality','assembly']:

        df['BASENAME'] = df.apply('{STRAIN}_{RUN}'.format_map, axis=1)

        df['ILLUMINA-R1'] = df['ILLUMINA_R1'].apply(lambda x:  os.path.abspath(x) if pd.notna(x) else x )
        df['ILLUMINA-R2'] = df['ILLUMINA_R2'].apply(lambda x:  os.path.abspath(x) if pd.notna(x) else x )
        df['NANOPORE'] =  df['NANOPORE'].apply(lambda x:  os.path.abspath(x) if pd.notna(x) else x )

        df['ILLUMINA_R1_RAW'] = df.apply("reads/raw/{BASENAME}_ILLUMINA-R1.fastq.gz".format_map, axis=1)
        df['ILLUMINA_R2_RAW'] = df.apply("reads/raw/{BASENAME}_ILLUMINA-R2.fastq.gz".format_map, axis=1)
        df['NANOPORE_RAW'] = df.apply("reads/raw/{BASENAME}_NANOPORE.fastq.gz".format_map, axis=1)

        df['QC_FILENAME_ILLUMINA_R1'] = df.apply("reads/qc/{BASENAME}_ILLUMINA-R1_fastqc.html".format_map,axis=1)
        df['QC_FILENAME_ILLUMINA_R2'] = df.apply("reads/qc/{BASENAME}_ILLUMINA-R2_fastqc.html".format_map,axis=1)
        df['QC_FILENAME_NANOPORE'] = df.apply("reads/qc/{BASENAME}_NANOPORE_nanostat.txt".format_map ,axis=1)

        df['ASSEMBLY_FILENAME'] =  df.apply("assembly/{STRAIN}_{RUN}.fasta".format_map, axis=1)


        QC = [list(df.dropna(subset=['ILLUMINA_R1'])['QC_FILENAME_ILLUMINA_R1']),
              list(df.dropna(subset=['ILLUMINA_R2'])['QC_FILENAME_ILLUMINA_R2']),
              list(df.dropna(subset=['NANOPORE'])['QC_FILENAME_NANOPORE']),
              ]


        config['DF'] = df

        #df = df.melt().dropna()

        print(df['ASSEMBLY_FILENAME'])
        return QC, df['ASSEMBLY_FILENAME']


    ##################################
    ##      contigs as input        ##
    ##################################

    elif mode in ['annotation']:

        df['BASENAME'] = df.apply('{GENUS}_{SPECIES}_{STRAIN}_{RUN}'.format_map, axis=1)
        config['DF'] = df
        ANNOTATIONS = ['abricate', 'rgi', 'mefinder', 'resfinder', 'virulencefinder', 'bacmetscan',
                       'islandpath', 'phispy', 'plasmidfinder', 'crisprcasfinder','keyword',
                       ]

        return [ f'{basename}.{program}.gff'
                 for basename in df.apply('annotation/special/{BASENAME}'.format_map,axis=1)
                 for program in ANNOTATIONS
                 ], \
               df.apply('annotation/special/mlst/{BASENAME}.mlst.txt'.format_map,axis=1),\
               df.apply('annotation/{BASENAME}.merged.gbk'.format_map, axis=1)

rule all:
    input:
        targets(pd.read_csv(config['csv']), config['mode'])

#concatenate and copy files into working directory
module ImportFiles:
    snakefile: r"rules/ImportFiles.smk"
    config: config
use rule * from ImportFiles as IN_*

# qc using nanostat and fastqc
module QualityControl:
    snakefile: r"rules/QualityControl.smk"
    config: config
use rule * from QualityControl as QC_*

module ErrorCorrection:
    snakefile: r"rules/ErrorCorrection.smk"
    config: config
use rule * from ErrorCorrection as EC_*

module Assembly:
    snakefile: r"rules/Assembly.smk"
    config: config
use rule * from Assembly as AS_*

module Annotation:
    snakefile: r"rules/Annotation.smk"
    config: config
use rule * from Annotation as AN_*