from subprocess import call
print("loading libs")


from os.path import join as pjoin
import os
import shutil

print("loading configs")

configfile : "repos/metasssnake/snakefiles/params.json"

config['general']['temp_dir'] = os.environ['SNIC_TMP']

print("loading samples")

samples = list(set(glob_wildcards("0000_raws/0100_reads/genomic2/{sample}_{lib}_R{dir}.fastq.gz").sample ))
print("We have ", len(samples), "samples")

coasses = glob_wildcards("9000_metadata/9100_samplesets/{coass}.txt").coass
print("We have ", len(coasses), "coasses")

print("loading snakefiles")


include : "read_processing.py"
include : "assembling.py"
include : "mapping.py"
include : "binning.py"
include : "mags.py"
include : "analyses.py"

print("loading shell")

shell.prefix("module load bioinfo-tools bbmap samtools BioPerl prokka  perl_modules; ")

localrules: all,clean_metabat,filter_assembly,hmmer_table, MAG_stats,full_kaiju_by_level,all_kaijus, mags2cogs

"""
DIR=`pwd`
QSCRIPT="sbatch -D $DIR -A snic2019-3-22 -t {cluster.time} -n {cluster.n} -p {cluster.part} -C {cluster.C} -M {cluster.M} -o {cluster.out} -e {cluster.err}"

snakemake -j 999 --cluster-config 8000_scripts/cluster.json --cluster "$QSCRIPT"  -s repos/metasssnake/snakefiles/main.py
"""


#temp = glob_wildcards("1000_processed_reads/{sample}/fwd.fastq.gz").sample

print("loading The Rule")

rule all : 
    input : expand("/crex/proj/uppstore2018116/moritz6/1000_processed_reads/{name}/assemblies/{assembler}/binning/metabat/full_taxonomy.tax", name = samples, assembler = ["megahit"]), expand("/crex/proj/uppstore2018116/moritz6/1500_coasses/{name}/assemblies/{assembler}/binning/metabat/abundance_tables/abundance_per_bin.csv", name = [a for a in coasses if a!="Loclat"  and a != "MJ-time" ], assembler = ["megahit"])



#)


