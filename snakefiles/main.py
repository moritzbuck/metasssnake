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

localrules: make_all,clean_metabat,filter_assembly,hmmer_table,make_all, MAG_stats, bbmap_index

"""
DIR=`pwd`
QSCRIPT="sbatch -D $DIR -A snic2019-3-22 -t {cluster.time} -n {cluster.n} -p {cluster.part} -C {cluster.C} -M {cluster.M} -o {cluster.out} -e {cluster.err}"

snakemake -j 999 --cluster-config 8000_scripts/cluster.json --cluster "$QSCRIPT"  -s repos/metasssnake/snakefiles/main.py
"""

def all_binnings(wildcards):
    samps = [f.split("_")[0] for f in os.listdir("0000_raws")]
    coass = [f.split(".")[0] for f in os.listdir("9000_metadata/9100_samplesets/")]

    proto = "{root}/{name}/assemblies/{assembler}/binning/metabat/normed_pfam_covs.csv"
    proto2 = "{root}/{name}/assemblies/{assembler}/binning/metabat/abundance_tables/abundance_per_bin.csv"
    singles = [proto.format(root = "/crex/proj/uppstore2018116/moritz6/1000_processed_reads", name = s, assembler = ass) for s in samps for ass in ['megahit']  if "WRT" in s]#, 'spades']]
    ningles = [proto.format(root = "/crex/proj/uppstore2018116/moritz6/1500_coasses", name = s, assembler = ass) for s in coass for ass in ['megahit'] ]
    ningles2 = [proto2.format(root = "/crex/proj/uppstore2018116/moritz6/1500_coasses", name = s, assembler = ass) for s in coass for ass in ['megahit'] ]
    return ningles2# +ningles #+ singles 

rule make_all:
    input : all_binnings
    output : "done"
    run :
        call("touch done", shell= True)


