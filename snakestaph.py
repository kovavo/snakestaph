#!/usr/bin/env python3
'''
Convience script to prepare and run snakestaph snakemake workflow
Author Vojtech Kovarovic kovarovic@mail.muni.cz
'''

import os
import argparse
import snakemake


#argument parsing

parser = argparse.ArgumentParser(description='Snakemake workflow to assembly, annotate and analyse staphylococcal genomes')

parser.add_argument('mode', help='Choose workflow', choices=['assembly', 'quality', 'annotation'])

parser.add_argument('-i','--input', help='Text file with columns') #TODO: define format
parser.add_argument('-o','--output', help='Set output directory', required=True)

#set location of conda envs and scratch dir, helpful in grid computing where more physical computers are used
parser.add_argument('-e', '--enviroments', help='Directory where the conda and singulaity enviroments will be stored', default=os.path.join(os.getcwd(),'snakestaph_env'))
parser.add_argument('-s', '--shadow', help='Conda shadow directory, in grid computing it should be set to scratch directory', default=os.path.join(os.getcwd(),'snakestaph_shadow'))

parser.add_argument('-t', '--threads', help='number of threads used by snakemake', default=8)
parser.add_argument('-u', '--unlock', help='Unlock snakemake directory (eg. after unexpected power off)', action='store_true', default=False)
parser.add_argument('-d', '--dryrun', action='store_true') #this is flag

args = parser.parse_args()

#set paths
snakefile = os.path.join(os.path.dirname(__file__), 'workflow/Snakefile.smk')
configfile = [os.path.join(os.path.dirname(__file__), 'config/config.yaml')]


#execute main snakemake
snakemake.snakemake(snakefile,
                    # basic arguments i like
                    keepgoing=True,
                    printshellcmds=True,

                    # switchs useful for debugging
                    unlock=args.unlock,
                    dryrun=args.dryrun,
                    #TODO: add overrides as string? Something like conda_cmd = ''
                    #verbose=True,
                    #debug_dag=True,
                    #printrulegraph=True,
                    #forceall=True,
                    #printdag=True,
                    #cleanup_metadata='annotation/special/mge/Staphylococcus_brunensis_NRL19-737.mge.gff',

                    # snakemake configuration
                    use_conda=True,
                    use_singularity=True,
                    shadow_prefix=args.shadow,
                    conda_prefix=os.path.abspath(args.enviroments),
                    singularity_prefix=os.path.abspath(args.enviroments),
                    conda_frontend='mamba',
                    workdir=args.output,
                    cores=int(args.threads),
                    configfiles=configfile,
                    config={'csv': os.path.abspath(args.input),
                            'mode': args.mode,
                            'database_prefix': args.enviroments, #for sequence db to be dowloaded for annotations
                            },
                    )