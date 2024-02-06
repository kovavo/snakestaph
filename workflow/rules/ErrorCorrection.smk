'''
pipeline for error correcting illumina, iontorrent and nanopore reads
multiple rules to encounter irregular combinations
allowed combinations are illumina+nanopore, iontorrent+nanopore, nanopore, illumina, iontorrent
'''
#prefer using short read reference if possible
ruleorder:  filtlong_short_read_reference > filtlong > fastp

rule fastp:
    input:
        left = lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['ILLUMINA_R1_RAW'],
        right = lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['ILLUMINA_R2_RAW'],
    output:
        left = 'reads/corrected/{filename}_ILLUMINA_R1.fastq.00.0_0.cor.fastq.gz',
        right = 'reads/corrected/{filename}_ILLUMINA_R2.fastq.00.0_0.cor.fastq.gz',
        single  = 'reads/corrected/{filename}_ILLUMINA_R_unpaired.00.0_0.cor.fastq.gz',

    threads: workflow.cores
    conda: '../envs/fastp.yaml'
    shadow: 'minimal'
    message:
        '\nError correction using spades fastp: {input} \n'
    shell:
        'fastp --thread {threads} -i {input.left} -I {input.right} -o {output.left} -O {output.right} --unpaired1 {output.single}  --unpaired2 {output.single}'

# filtering nanopore reads
rule filtlong:
    input:
        lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['NANOPORE_RAW'],
    output:
        'reads/corrected/{filename}_NANOPORE.fastq.gz'
    threads: workflow.cores
    conda: '../envs/filtlong.yml'
    params:
        min_length = 100,
        keep_percent = 99,
        target_bases =  1000, #count in mb
    message:
        '\nFilter reads based on quality and length: {input} \n'
    shell:
        'filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} --target_bases {params.target_bases} {input} | gzip > {output}'

#use paired end from illumina for reference
rule filtlong_short_read_reference:
    input:
        long = lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['NANOPORE_RAW'],
        left = lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['ILLUMINA_R1_RAW'],
        right = lambda wildcards: config['DF'].loc[config['DF']['READ_BASENAME'] == wildcards.filename]['ILLUMINA_R2_RAW'],
    output:
        'reads/corrected/{filename}_NANOPORE-SE.fastq.gz'
    threads: workflow.cores
    conda: '../envs/filtlong.yml'
    params:
        min_length = 100,
        keep_percent = 99,
        target_bases =  1000, #count in mb
    message:
        '\nFilter reads based on quality and length: {input} \n'
    shell:
        'filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} --target_bases {params.target_bases} '
        '{input.long} -1 {input.left} -2 {input.right} | gzip > {output}'
