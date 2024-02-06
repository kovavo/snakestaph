'''
Annotate bacterial genome
'''

#TODO: add some BGC calling pipeline
#TODO: add PGAP pipeline - metacentrum only too large
#TODO: convert this code to be run on metacentrum!!!

import glob


#setup databases
ruleorder: islandpath > phispy > keyword > convert_to_gff > merge

###################################
##       general annotation      ##
###################################

rule prokka:
    input:
        rules.IN_fasta_assemblies.output
    output:
        # maybe not necessary to write down all output files, since the shadow is not used
        errors = "annotation/prokka/{sample}.prokka.err",
        protein= "annotation/prokka/{sample}.prokka.faa",
        fasta= "annotation/prokka/{sample}.prokka.fna",
        ffn = "annotation/prokka/{sample}.prokka.ffn",
        fsa = "annotation/prokka/{sample}.prokka.fsa",
        genbank = "annotation/prokka/{sample}.prokka.gbk",
        gff = "annotation/prokka/{sample}.prokka.gff",
        log = "annotation/prokka/{sample}.prokka.log",
        sequin = "annotation/prokka/{sample}.prokka.sqn",
        tbl = "annotation/prokka/{sample}.prokka.tbl",
        table = "annotation/prokka/{sample}.prokka.tsv",
        stats = "annotation/prokka/{sample}.prokka.txt",
    params:
        kingdom = 'Bacteria',
        genus= lambda wildcards: config['DF'].loc[config['DF']['BASENAME'] == wildcards.sample]['GENUS'][0],
        strain = lambda wildcards: config['DF'].loc[config['DF']['BASENAME'] == wildcards.sample]['STRAIN'][0],
        species = lambda wildcards: config['DF'].loc[config['DF']['BASENAME'] == wildcards.sample]['SPECIES'][0],
    conda: '../envs/prokka.yaml'
    threads: workflow.cores
    priority: 2 #run first
    shell:
         "gunzip {input}; "
         "prokka --force --usegenus --cpus {threads} --genus {params.genus} --species {params.species} --kingdom {params.kingdom} "
         "--strain {params.strain} --prefix {wildcards.sample}.prokka --outdir annotation/prokka/ annotation/assembly/{wildcards.sample}.fasta --locustag cluster --increment 1"

###################################
##   virulence and resistance    ##
###################################

rule abricate:
    input:
        rules.prokka.output.fasta,
    output:
        "annotation/special/abricate/{sample}.abricate"
    threads: workflow.cores
    conda: '../envs/abricate.yaml'
    #shadow: 'minimal'
    params:
        databases= ['vfdb', 'card', 'argannot', 'plasmidfinder', 'resfinder', 'ncbi', 'megares']
    shell:
         "for DATABASE in {params.databases}; do "
         "abricate --threads {threads} --db $DATABASE {input} > annotation/special/abricate/{wildcards.sample}.abricate-${{DATABASE}}.tab; "
         "done ;"
         "awk FNR-1 annotation/special/abricate/{wildcards.sample}*.tab > {output}"

rule virulencefinder:
    input:
        fasta=rules.prokka.output.fasta,
    output:
        'annotation/special/virulencefinder/{sample}.virulencefinder'
    threads: workflow.cores
    conda: '../envs/virulencefinder.yaml'
    params:
        db_path=config['database_prefix']
    #shadow: 'minimal'
    message:
        '\nRun virulencefinder for \n'
    shell:
        'cd {params.db_path};'
        'download-virulence-db.sh;'
        'cd -;'
        
        'virulencefinder.py --databasePath {params.db_path}/virulencefinder_db '
        '-i {input.fasta} -o annotation/special/virulencefinder;'
        
        'mv annotation/special/virulencefinder/data.json {output}'

rule rgi:
    input:
        rules.prokka.output.fasta
    output:
        'annotation/special/rgi/{sample}.rgi',
    threads: workflow.cores
    conda: '../envs/rgi.yaml'
    #shadow: 'minimal'
    message:
        '\nRun rgi for \n'
    shell:
        'rgi main --input_type contig --clean --num_threads {threads} '
        '--input_sequence {input} --output_file annotation/special/rgi/{wildcards.sample}.rgi;'
        'mv annotation/special/rgi/{wildcards.sample}.rgi.txt {output}'

rule resfinder:
    input:
        fasta = rules.prokka.output.fasta,
    output:
        'annotation/special/resfinder/{sample}.resfinder'
    threads: workflow.cores
    conda: '../envs/resfinder.yaml'
    params:
        db_path=config['database_prefix']
    #shadow: 'minimal'
    message:
        '\nRun resfinder for \n'
    shell:
        'RESFINDER_DB={params.db_path};'
        'download-db.sh;'
        'python -m resfinder -s other --acquired -db_res {params.db_path}/db_resfinder -ifa {input.fasta} -o annotation/special/resfinder;'
        'mv annotation/special/resfinder/ResFinder_results_tab.txt {output}'

rule bacmetscan:
    input:
        protein = rules.prokka.output.protein,
    output:
        'annotation/special/bacmetscan/{sample}.bacmetscan'
    threads: workflow.cores
    conda: '../envs/bacmetscan.yaml'
    #shadow: 'minimal'
    message:
        '\nRun bacmetscan for \n'
    params:
        db_path=config['database_prefix'],
        database='EXP'
    shell:
        'if test -d {params.db_path}/db_bacmet2;'
        '  then echo "BacMet2 database already downloaded";'
        '  else mkdir {params.db_path}/db_bacmet2;'
        '  cd {params.db_path}/db_bacmet2;'
        '  wget http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta ;'
        '  wget http://bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz; '
        '  wget http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt;'
        '  mv BacMet2_EXP.753.mapping.txt BacMet2_EXP_database.mapping.txt;' #rename mapping file so bacmetscan recognize it
        '  wget http://bacmet.biomedicine.gu.se/download/BacMet-Scan;'
        '  cd;'
        'fi;'
        
        'perl {params.db_path}/db_bacmet2/BacMet-Scan -v -columns all -protein -cpu {threads} '
        '-d {params.db_path}/db_bacmet2/BacMet2_{params.database}_database -i {input.protein} -o annotation/special/bacmetscan;'
        'mv annotation/special/bacmetscan.table {output}'


###################################
##   virulence and resistance    ##
###################################


#mobile genetic elements
rule MobileGenomeElements:
    input:
        rules.prokka.output.fasta
    output:
        "annotation/special/mefinder/{sample}.mefinder"
    threads: workflow.cores
    container: 'docker://mkhj/mobile_element_finder:latest'
    message:
        '\nRun MobileGenomeElements for \n'
    shadow: 'minimal'
    shell:
        'mefinder find --gff --threads {threads} --contig {input} annotation/special/mefinder/{wildcards.sample}.mefinder;'
        'cp annotation/special/mefinder/{wildcards.sample}.mefinder.gff annotation/special/{wildcards.sample}.mefinder.gff;'
        'mv annotation/special/mefinder/{wildcards.sample}.mefinder.gff {output}'

rule islandpath:
    input:
        rules.prokka.output.genbank
    output:
        'annotation/special/{sample}.islandpath.gff'
    threads: workflow.cores
    conda: '../envs/islandpath.yaml'
    #shadow: 'minimal'
    message:
        '\nRun islandpath for \n'
    shell:
        'islandpath {input} annotation/special/islandpath/{wildcards.sample}.islandpath;'
        'cp annotation/special/islandpath/{wildcards.sample}.islandpath {output}'

rule phispy:
    input:
        rules.prokka.output.genbank,
    output:
        'annotation/special/{sample}.phispy.gff'
    threads: workflow.cores
    conda: '../envs/phispy.yaml'
    #shadow: 'minimal'
    message:
        '\nRun phispy for \n'
    shell:
        'PhiSpy.py {input} --output_choice 512 --nonprophage_genegaps 20 --threads {threads} -p {wildcards.sample}.phispy -o annotation/special/phispy;'
        "cp annotation/special/phispy/{wildcards.sample}.phispy_prophage.gff3 {output};"

rule plasmidfinder:
    input:
        fasta = rules.prokka.output.fasta,
    output:
        'annotation/special/plasmidfinder/{sample}.plasmidfinder'
    threads: workflow.cores
    conda: '../envs/plasmidfinder.yaml'
    params:
        db_path=config['database_prefix']
    #shadow: 'minimal'
    message:
        '\nRun plasmidfinder for \n'
    shell:
        'PLASMID_DB={params.db_path}/plasmidfinder_db;'
        'mkdir $PLASMID_DB;'
        'download-db.sh;'
        'plasmidfinder.py  -p {params.db_path}/plasmidfinder_db -i {input.fasta} -o annotation/special/plasmidfinder ;'
        'mv annotation/special/plasmidfinder/data.json {output}'

rule crisprcasfinder:
    input:
        rules.prokka.output.fasta
    output:
        'annotation/special/crisprcasfinder/{sample}.crisprcasfinder'
    threads: workflow.cores
    container: 'docker://unlhcc/crisprcasfinder:latest'
    #shadow: 'minimal'
    shell:
        'CRISPRCasFinder.pl -cas -in {input} -out annotation/special/crisprcasfinder/{wildcards.sample} -so /opt/CRISPRCasFinder/sel392v2.so;'
        'mv annotation/special/crisprcasfinder/{wildcards.sample}/result.json {output}'



####################################
##         other features         ##
####################################


rule keyword:
    input:
         rules.prokka.output.genbank
    output:
         "annotation/special/{sample}.keyword.gff"
    conda: '../envs/biopython.yaml'
    threads: workflow.cores
    params:
        feature = ['CDS'],
        qualifier = ['product'],
        values = ['phage','sapi','integrase','mobile','cassette','crispr','rlmh','resolvase','recombinase',
                   'transposase','capsular','modification','restriction','resistance','dnaA',
                  ]
    script:
         "../scripts/genbank_keywords.py"

#strain typing
rule mlst:
    input:
         rules.prokka.output.fasta
    output:
         "annotation/special/mlst/{sample}.mlst.txt"
    conda: '../envs/mlst.yaml'

    threads: workflow.cores
    shell:
         "mlst --threads {threads} {input} > {output}"

# #parse outputs
# #TODO: directory structure and mapping for parsing
# #TODO: change this so the functions are called inside of rules
#
checkpoint convert_to_gff:
    input:
        genbank = rules.prokka.output.genbank,
        outputs = "annotation/special/{program}/{sample}.{program}",
    output: 'annotation/special/{sample}.{program}.gff'
    threads: workflow.cores
    conda: '../envs/biopython.yaml'
    message:
        '\nRun parse for {wildcards.sample} \n'
    script:
        '../scripts/{wildcards.program}.py'

#TODO: change it to conditional checkpoint

rule merge:
    input:
        genbank = rules.prokka.output.genbank,
        gff = glob.glob('annotation/special/*.gff')
    output:
        "annotation/{sample}.merged.gbk"
    threads: workflow.cores
    conda: '../envs/biopython.yaml'
    message:
        '\nRun merge for \n'
    priority: -1 #run as last
    script:
        '../scripts/merge_genbank.py'

