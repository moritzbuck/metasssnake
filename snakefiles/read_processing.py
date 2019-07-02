
from subprocess import call
from os.path import join as pjoin
import os
import shutil

temp_dir = os.environ['SNIC_TMP']


def find_libs(wildcards):
    raw_path = config['general']['raw_folder']
    libs = [f.split("_")[1] for f in os.listdir(raw_path) if f.startswith(wildcards.sample + "_") and "_R1" in f]
    return ["1000_processed_reads/{sample}/reads/trimmomatic/{lib}/fwd_paired.fastq.gz".format(lib = l, sample = wildcards.sample) for l in libs]

# def find_fastq(wildcards):
#     path = '0000_raws/0100_reads/genomic/'
#     result = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*.fastq.gz')) if "/" + wildcards.sample +"_" in y]
#     assert len(result) == 2, print(result)
#     return result
#
#
# def all_dones(wildcards):
#     path = '0000_raws/0100_reads/genomic/'
#     result = [pjoin("1000_processed_reads/","_".join(s.split("_")[:-4]),"done") for s in os.listdir(path) if s.endswith(".fastq.gz")  ]
#     return list(set(result))
#
#
# rule fastqc:
#     input :  find_fastq
#     output : "1000_processed_reads/{sample}/reads/fastqc"
#     threads : THREADS
#     shell:
#         """
#         fastqc -t {threads} -o {output} {input}
#         """


rule trimmomatic:
    """ QCing and cleaning reads """
    params : java_cmd = config["read_processing"]['trimmomatic']['java_cmd'],
             jar_file = config["read_processing"]['trimmomatic']['jar_file'],
             mem = config["read_processing"]['trimmomatic']['java_vm_mem'],
             options = config["read_processing"]['trimmomatic']['options'],
             processing_options = config["read_processing"]['trimmomatic']['processing_options'],
             temp_folder = temp_dir
    input :  fwd = pjoin(config['general']['raw_folder'], "{sample}_{lib}_R1.fastq.gz"),
             rev = pjoin(config['general']['raw_folder'],"{sample}_{lib}_R2.fastq.gz")
    threads : 20
    output : read1 = "1000_processed_reads/{sample}/reads/trimmomatic/{lib}/fwd_paired.fastq.gz",
             read2 = "1000_processed_reads/{sample}/reads/trimmomatic/{lib}/rev_paired.fastq.gz",
             read1U = "1000_processed_reads/{sample}/reads/trimmomatic/{lib}/fwd_unpaired.fastq.gz",
             read2U = "1000_processed_reads/{sample}/reads/trimmomatic/{lib}/rev_unpaired.fastq.gz",
    log : "1000_processed_reads/{sample}/reads/trimmomatic/{lib}/log"
    shell:
        """
        mkdir -p {params.temp_folder}/{wildcards.sample}_{wildcards.lib}/
        unpigz -c -p {threads} {input.fwd}  >  {params.temp_folder}/{wildcards.sample}_{wildcards.lib}/temp_R1.fastq
        unpigz -c -p {threads} {input.rev} >  {params.temp_folder}/{wildcards.sample}_{wildcards.lib}/temp_R2.fastq

        {params.java_cmd} -Xmx{params.mem} -Djava.io.tmpdir={params.temp_folder} -jar {params.jar_file} PE {params.options} {params.temp_folder}/{wildcards.sample}_{wildcards.lib}/temp_R1.fastq {params.temp_folder}/{wildcards.sample}_{wildcards.lib}/temp_R2.fastq -threads {threads} {output.read1} {output.read1U} {output.read2}  {output.read2U} {params.processing_options} 2> {log}
        """

rule merge_libs:
    input : find_libs
    output : read1 = "1000_processed_reads/{sample}/reads/fwd.fastq.gz",
             read2 = "1000_processed_reads/{sample}/reads/rev.fastq.gz",
             unpaired = "1000_processed_reads/{sample}/reads/unp.fastq.gz"
    params : temp_folder = temp_dir
    threads : 20
    run :
        out_fold = pjoin(params.temp_folder, wildcards.sample)
        dirs = [os.path.dirname(l) for l in input]
        if not os.path.exists(out_fold):
            os.makedirs(out_fold)
        unzip_cmd = "unpigz -c -p {threads} {files} > {temp_fold}"
        call(unzip_cmd.format(threads = threads, files = " ".join([pjoin(d, "fwd_paired.fastq.gz" ) for d in dirs]), temp_fold = pjoin(out_fold, "fwd.fastq")), shell = True)
        call(unzip_cmd.format(threads = threads, files = " ".join([pjoin(d, "rev_paired.fastq.gz" ) for d in dirs]), temp_fold = pjoin(out_fold, "rev.fastq")), shell = True)
        call(unzip_cmd.format(threads = threads, files = " ".join([pjoin(d, "fwd_unpaired.fastq.gz" ) for d in dirs] + [pjoin(d, "rev_unpaired.fastq.gz" ) for d in dirs]), temp_fold = pjoin(out_fold, "unp.fastq")), shell = True)
        call("pigz -p {threads} {temp_fold}/*.fastq".format(temp_fold = out_fold, threads = threads), shell = True)
        call("mv {temp_fold}/*.gz {outdir}/".format(temp_fold = out_fold, outdir = os.path.dirname(output.read1), threads = threads), shell = True)


