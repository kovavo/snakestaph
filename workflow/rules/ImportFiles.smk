'''
Module for importing source reads and assemblies
'''

rule fastq_reads:
    input:
        lambda wildcards: config['DF'].loc[config['DF']['BASENAME'] == wildcards.filename][wildcards.platform],
    output:
       "reads/raw/{filename}_{platform}.fastq.gz"
    message: "Import raw reads {output} from {input}"
    threads: 1

    script:
        '../scripts/import_file.py'


rule fasta_assemblies:
    input:
        lambda wildcards: config['DF'].loc[config['DF']['BASENAME'] == wildcards.sample]['PATH']
    threads: 1
    output:
        'annotation/assembly/{sample}.fasta.gz'
    message: "Import {output}"
    script:
        '../scripts/import_file.py'
