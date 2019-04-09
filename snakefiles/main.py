from subprocess import call
from os.path import join as pjoin
import os
import shutil

configfile : "repos/metasssnake/snakefiles/params.json"

config['general']['temp_dir'] = os.environ['SNIC_TMP']

include : "read_processing.py"
include : "assembling.py"
include : "mapping.py"
include : "binning.py"
include : "mags.py"

shell.prefix("module load bioinfo-tools bbmap samtools BioPerl prokka  perl_modules; ")


def all_binnings(wildcards):
    samps = [f.split("_")[0] for f in os.listdir("0000_raws")]
    coass = [f.split(".")[0] for f in os.listdir("9000_metadata/9100_samplesets/") if f != "FL.txt" and f != "WRT"]

    proto = "{root}/{name}/assemblies/{assembler}/binning/metabat/normed_pfam_covs.csv"

    singles = [proto.format(root = "1000_processed_reads", name = s, assembler = ass) for s in samps for ass in ['megahit']  if "WRT" in s]#, 'spades']]
    ningles = [proto.format(root = "1500_coasses", name = s, assembler = ass) for s in coass for ass in ['megahit']]#, 'spades']]
    return ningles #+ singles 

rule make_all:
    input : all_binnings
    output : "done"
    run :
        call("touch done", shell= True)


