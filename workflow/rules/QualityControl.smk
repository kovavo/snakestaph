'''
Quality control and statistics for sequencing reads
'''

#QC and basic stats for both illumina and iontorrent reads
#"0.74.0/bio/fastqc

rule FastQC:
    input:
        #raw_reads = "reads/raw/{filename}_{platform}.fastq.gz"
        raw_reads = rules.IN_fastq_reads.output
    output:
        fastqc_html = "reads/qc/{filename}_{platform}_fastqc.html"

    threads: workflow.cores*config['paralel']*0.5
    priority: 1
    conda: "../envs/fastqc.yaml"
    shadow: "shallow"


    message:
        '\nQuality control with FastQC on following reads: {input}\n'
        '{threads} threads avaiable'
    shell:
        'fastqc -t {threads} -o reads/qc {input.raw_reads}'

#basic stats for nanopore reads
rule nanostat:
    input:
        raw_reads= rules.IN_fastq_reads.output

    output:
        "reads/qc/{filename}_{platform}_nanostat.txt"
    threads: workflow.cores*config['paralel']*0.5
    conda: '../envs/nanostat.yaml'
    shadow: "shallow"
    message:
        '\nComputing read statistics using NanoStat: {input} \n'
    shell:
        'NanoStat -t {threads} --tsv --fastq {input} > {output}'
