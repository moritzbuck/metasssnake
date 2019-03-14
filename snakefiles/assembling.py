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

rule filter_assembly:
    params : cutoff = config['assembling']['filter_cutoff']
    input : assembly = "{path}/{sample}/assemblies/{assembler}/assembly.fna",
    output : filtered_assembly = "{path}/{sample}/assemblies/{assembler}/filtered_assembly.fna",
    threads : 1
    run :
        from Bio import SeqIO
        out_contigs = []
        count = 1
        for s in SeqIO.parse(input.assembly, "fasta"):
            if len(s.seq) >= params.cutoff:
                s.id = "{type}_{assembler}_{sample}_{num:07}".format(type = "coass" if "1500_" in wildcards.path else "single", sample = wildcards.sample, assembler = wildcards.assembler, num = count)
                s.description = ""
                count += 1
                out_contigs += [s]
        SeqIO.write(out_contigs, output.filtered_assembly, "fasta")
                

rule assemble:
    params : temp_folder = pjoin(config['general']['temp_dir'], "{sample}", "{assembler}"),
    input : unpack(get_libs),
    output : assembly = "{path}/{sample}/assemblies/{assembler}/assembly.fna",
    threads : 16,
    run :
        os.makedirs(params.temp_folder, exist_ok = True)
        unzip_cmd = "unpigz -c -p {threads} {files} > {temp_fold}"
        fwd = pjoin(params.temp_folder, "fwd.fastq")
        rev = pjoin(params.temp_folder, "rev.fastq")
        unp = pjoin(params.temp_folder, "unp.fastq")

        if os.path.exists(pjoin(os.path.dirname(output.assembly),"reads", "fwd.bbnorm.fastq.gz")) :
            call(unzip_cmd.format(threads = threads, files = pjoin(os.path.dirname(output.assembly),"reads", "fwd.bbnorm.fastq.gz"), temp_fold = fwd), shell = True)
            call(unzip_cmd.format(threads = threads, files = pjoin(os.path.dirname(output.assembly),"reads", "rev.bbnorm.fastq.gz"), temp_fold = rev), shell = True)
#            call("touch {files}". format(files = unp), shell = True)
        else :
            call(unzip_cmd.format(threads = threads, files = " ".join(input.fwd) if type(input.fwd) == list else input.fwd , temp_fold = fwd), shell = True)
            call(unzip_cmd.format(threads = threads, files = " ".join(input.rev) if type(input.rev) == list else input.rev , temp_fold = rev), shell = True)
            call(unzip_cmd.format(threads = threads, files = " ".join(input.unp) if type(input.unp) == list else input.unp , temp_fold = unp), shell = True)

        if wildcards.assembler == "megahit":
            if os.path.exists(unp):
                call("megahit -m 0.90 -1 {fwd} -2 {rev} -r {unp} -t {threads} -o {outfold} --out-prefix megahit".format(fwd = fwd, rev = rev, unp = unp, threads = threads, outfold = pjoin(params.temp_folder, "data")), shell = True)
            else :
                call("megahit -m 0.90 -1 {fwd} -2 {rev} -t {threads} -o {outfold} --out-prefix megahit".format(fwd = fwd, rev = rev, threads = threads, outfold = pjoin(params.temp_folder, "data")), shell = True)
            shutil.rmtree(pjoin(params.temp_folder, "data", "intermediate_contigs"))
            shutil.move(pjoin(params.temp_folder, "data"), pjoin(os.path.dirname(output.assembly),"data"))
            os.symlink(pjoin(os.getcwd(),pjoin(os.path.dirname(output.assembly),"data"), "megahit.contigs.fa"), output.assembly)

        elif wildcards.assembler == "spades":
            if os.path.exists(unp):
                call("spades.py --meta  -1 {fwd} -2 {rev} -s {unp} -t {threads} -o {outfold}".format(fwd = fwd, rev = rev, unp = unp, threads = threads, outfold = params.temp_folder), shell = True)
            else :
                call("spades.py --meta  -1 {fwd} -2 {rev} -t {threads} -o {outfold}".format(fwd = fwd, rev = rev, threads = threads, outfold = params.temp_folder), shell = True)
#            shutil.rmtree(pjoin(params.temp_folder, intermediate_contigs))
            shutil.rmtree(pjoin(params.temp_folder, "K21"))
            shutil.rmtree(pjoin(params.temp_folder, "K33"))
            shutil.rmtree(pjoin(params.temp_folder, "K55"))
            shutil.move(params.temp_folder, pjoin(os.path.dirname(output.assembly),"data") )
            os.symlink(pjoin(os.getcwd(),pjoin(os.path.dirname(output.assembly),"data"), "scaffolds.fasta"), output.assembly)
        else :
            print("Not an accepted assembler")
            return "broken"
