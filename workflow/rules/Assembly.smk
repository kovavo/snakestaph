'''
pipeline for illumina and nanopore hybird assembly
may need some parameter tinkering
'''

#prefer hybrid assembly if possible
ruleorder: unicycler_hybrid > unicycler_short_reads > unicycler_long_reads

rule unicycler_hybrid:
    input:
        long=rules.EC_filtlong_short_read_reference.output,
        left=rules.EC_fastp.output.left,
        right=rules.EC_fastp.output.right,
        single=rules.EC_fastp.output.single,
    output:
        graph='assembly/unicycler/{filename}/assembly.gfa',
        contigs='assembly/{filename}.fasta',
    threads: workflow.cores
    conda: '../envs/unicycler.yaml'
    #shadow: 'minimal'
    params:
        params='--mode bold --min_fasta_length 0 --verbosity 2 --start_gene_id 10 --start_gene_cov 10 --kmers 33,55,77,99,115,125,127',
    message:
        '\nLong read only assembly: {input} \n'
    shell:
        'unicycler -t {threads} -l {input.long} -1 {input.left} -2 {input.right} -s {input.single} '
        '-o assembly/unicycler/{wildcards.filename} {params.params} ;'

        #rename output contigs and move them to single directory
        'cp assembly/unicycler/{wildcards.filename}/assembly.fasta '
        'assembly/{wildcards.filename}.fasta'


rule unicycler_short_reads:
    input:
        left=rules.EC_fastp.output.left,
        right=rules.EC_fastp.output.right,
        single=rules.EC_fastp.output.single,

    output:
        graph='assembly/unicycler/{filename}/assembly.gfa',
        contigs='assembly/{filename}.fasta',
    threads: workflow.cores
    conda: '../envs/unicycler.yaml'
    #shadow: 'minimal'
    params:
        params='--mode bold --min_fasta_length 0 --verbosity 2 --start_gene_id 10 --start_gene_cov 10 --kmers 33,55,77,99,115,125,127',

    message:
        '\nLong read only assembly: {input} \n'
    shell:
        'unicycler -t {threads} -1 {input.left} -2 {input.right} -s {input.single} '
        '-o assembly/unicycler/{wildcards.filename}  {params.params}  ;'
        
         #rename output contigs and move them to single directory
        'cp assembly/unicycler/{wildcards.filename}/assembly.fasta '
        'assembly/{wildcards.filename}.fasta'

rule unicycler_long_reads:
    input:
        rules.EC_filtlong.output,
    output:
        graph='assembly/unicycler/{filename}/assembly.gfa',
        contigs='assembly/{filename}.fasta'
    threads: workflow.cores
    conda: '../envs/unicycler.yaml'
    params:
        params = '--mode bold --min_fasta_length 0 --verbosity 2 --start_gene_id 10 --start_gene_cov 10 --kmers 33,55,77,99,115,125,127',
    message:
        '\nLong read only assembly: {input} \n'
    shell:
        'unicycler -t {threads} -l {input} -o assembly/unicycler/{wildcards.filename} {params.params} ;'

        #rename output contigs and move them to single directory
        'cp assembly/unicycler/{wildcards.filename}/assembly.fasta '
        'assembly/{wildcards.filename}.fasta'
