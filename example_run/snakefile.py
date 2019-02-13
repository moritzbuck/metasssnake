
from subprocess import call
from os.path import join as pjoin
import os
import shutil

configfile : "/home/moritz/people/0023_anoxicencyclo/8000_scripts/snakemake_config.json"

config['general']['temp_dir'] = os.environ['SNIC_TMP']

include : "/home/moritz/repos/moritz/metasssnake/snakefiles/read_processing.py"
include : "/home/moritz/repos/moritz/metasssnake/snakefiles/assembling.py"
include : "/home/moritz/repos/moritz/metasssnake/snakefiles/mapping.py"
include : "/home/moritz/repos/moritz/metasssnake/snakefiles/binning.py"
include : "/home/moritz/repos/moritz/metasssnake/snakefiles/mags.py"

shell.prefix("module load bioinfo-tools bbmap samtools perl_modules BioPerl prokka; ")


def all_binnings(wildcards):
    samps = set([f.split("_")[0] for f in os.listdir("0000_raws/0100_reads/genomic2/") if "LaPlata" in f])
    coass = [f.split(".")[0] for f in os.listdir("9000_metadata/9100_samplesets/")]

    proto = "{root}/{name}/assemblies/{assembler}/binning/metabat/magstats.csv"

    singles = [proto.format(root = "1000_processed_reads", name = s, assembler = ass) for s in samps for ass in ['megahit']]#, 'spades']]
#    ningles = [proto.format(root = "1500_coasses", name = s, assembler = ass) for s in coass for ass in ['megahit']]#, 'spades']]
    return singles 

rule make_all:
    input : all_binnings
    output : "done"
    run :
        call("touch done", shell= True)


