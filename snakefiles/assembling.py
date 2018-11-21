from subprocess import call
from os.path import join as pjoin
import os
import shutil



def get_libs(wildcards):
    sets_dir = pjoin(config['general']['home'] , "9000_metadata/9100_samplesets/")
    sample_sets = [ f[:-4] for f in os.listdir(sets_dir) ]
    if wildcards.sample in sample_sets:
        with open(pjoin(sets_dir, wildcards.sample +".txt")) as handle:
            samples = [l.rstrip() for l in handle]
        output = {
            'fwd' : ["1000_processed_reads/{sample}/reads/fwd.fastq.gz".format(sample = s) for s in samples],
            'rev' : ["1000_processed_reads/{sample}/reads/rev.fastq.gz".format(sample = s) for s in samples],
            'unp' : ["1000_processed_reads/{sample}/reads/unp.fastq.gz".format(sample = s) for s in samples]}
        return output
    else : 
        output = {
            'fwd' : "1000_processed_reads/{sample}/reads/fwd.fastq.gz".format(sample = wildcards.sample),
            'rev' : "1000_processed_reads/{sample}/reads/rev.fastq.gz".format(sample = wildcards.sample),
            'unp' : "1000_processed_reads/{sample}/reads/unp.fastq.gz".format(sample = wildcards.sample)}
        return output


rule assemble:
    params : temp_folder = pjoin(config['general']['temp_dir'], "{sample}", "{assembler}")
    input : unpack(get_libs)
    output : assembly = "{path}/{sample}/assemblies/{assembler}/assembly.fna",
             folder = "{path}/{sample}/assemblies/{assembler}/data"
    threads : 20
    run :
        os.makedirs(params.temp_folder, exist_ok = True)
        unzip_cmd = "unpigz -c -p {threads} {files} > {temp_fold}"
        fwd = pjoin(params.temp_folder, "fwd.fastq")
        rev = pjoin(params.temp_folder, "rev.fastq")
        unp = pjoin(params.temp_folder, "unp.fastq")

        call(unzip_cmd.format(threads = threads, files = " ".join(input.fwd) if type(input.fwd) == list else input.fwd , temp_fold = fwd), shell = True)
        call(unzip_cmd.format(threads = threads, files = " ".join(input.rev) if type(input.rev) == list else input.rev , temp_fold = rev), shell = True)
        call(unzip_cmd.format(threads = threads, files = " ".join(input.unp) if type(input.unp) == list else input.unp , temp_fold = unp), shell = True)

        if wildcards.assembler == "megahit":
            call("megahit --continue -1 {fwd} -2 {rev} -r {unp} -t {threads} -o {outfold} --out-prefix megahit".format(fwd = fwd, rev = rev, unp = unp, threads = threads, outfold = pjoin(params.temp_folder, "data")), shell = True)
            shutil.rmtree(pjoin(params.temp_folder, "data", intermediate_contigs))
            shutil.move(pjoin(params.temp_folder, "data"), output.folder)
            os.symlink(pjoin(output.folder, "megahit.contigs.fa"), output.assembly)

        elif wildcards.assembler == "spades":
            call("spades.py --meta  -1 {fwd} -2 {rev} -s {unp} -t {threads} -o {outfold}".format(fwd = fwd, rev = rev, unp = unp, threads = threads, outfold = params.temp_folder), shell = True)
#            shutil.rmtree(pjoin(params.temp_folder, intermediate_contigs))
            shutil.move(params.temp_folder, output.folder)
            os.symlink(pjoin(output.folder, "scaffolds.fasta"), output.assembly)
        else :
            print("Not an accepted assembler") 
            return "broken"

