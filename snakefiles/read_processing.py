
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


rule mash:
    params : kmer = config["read_processing"]['mash']['kmer'],
             hashes = config["read_processing"]['mash']['hashes']
    input : "1000_processed_reads/{sample}/reads/fwd.fastq.gz","1000_processed_reads/{sample}/reads/rev_2P.fastq.gz"
    output : "1000_processed_reads/{sample}/reads/mash/{sample}.msh"
    log : "1000_processed_reads/{sample}/reads/mash/{sample}.log"
    threads : 1
    shell :
        "mash sketch -r -p {threads} -k {params.kmer} -s {params.hashes} -o $(echo '{output}' | sed -e 's/.msh//') {input} > {log}"

rule kaiju: 
    params : db_path = config["read_processing"]['kaiju']['db_path'],
    input : fwd = "1000_processed_reads/{sample}/reads/fwd.fastq.gz",  rev = "1000_processed_reads/{sample}/reads/rev.fastq.gz",
    output : "1000_processed_reads/{sample}/reads/kaiju/{db}/{sample}.kaiju"
    log : "1000_processed_reads/{sample}/reads/kaiju/{db}/kaiju.log"
    threads : 20
    shell : """
    kaiju -t {params.db_path}/{wildcards.db}/nodes.dmp -f {params.db_path}/{wildcards.db}/*.fmi   -i {input.fwd} -j {input.rev} -o {output} -z {threads} > {log}
    """

rule full_kaiju_by_level:
    params : db_path = config["read_processing"]['kaiju']['db_path'],
    input : expand("1000_processed_reads/{sample}/reads/kaiju/{{db}}/{sample}.kaiju", sample = samples)
    output : "4000_library_analysis/kaiju/{db}/tables/{level}.tsv"
    shell : """
kaiju2table -p -t {params.db_path}/{wildcards.db}/nodes.dmp -n {params.db_path}/{wildcards.db}/names.dmp -r {wildcards.level} -o {output} {input} 
"""

rule all_kaijus:
    input : expand("4000_library_analysis/kaiju/{db}/tables/{level}.tsv", db = ['nr_euk_mycocosm'], level = ['phylum', 'class', 'order', 'family', 'genus', 'species'])
    output : "4000_library_analysis/kaiju/kaiju.ed"
    shell: """
touch {output}
"""
#
# rule matam:
#     params : db_path = config["read_processing"]['matam']['db_path'],
#              temp_dir = config['temp_dir'],
#              max_mem = config["read_processing"]['matam']['max_mem']
#     threads : config['threads']
#     input :  fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
#              rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
#     output : dir =  "1000_processed_reads/{sample}/reads/matam/",
#              fasta =  "1000_processed_reads/{sample}/reads/matam/final_assembly.fa",
#     shell : """
#     paste <(unpigz  -p {threads} -c {input.fwd} | paste  - - - -) <(unpigz {threads} -c {input.rev} | paste - - - -) | tr $"\t" $"\n" >  {params.temp_dir}/matam_{wildcards.sample}.fastq
#     matam_assembly.py -d {params.db_path} -i {params.temp_dir}/matam_{wildcards.sample}.fastq --cpu {threads} --max_memory {params.max_mem} -v --perform_taxonomic_assignment -o {output.dir}
#     """
#
#
# rule all_kaiju:
#     input : "1000_processed_reads/done"
#     output : "1000_processed_reads/kaiju_table.tsv"
#     run :
#         from ete3 import NCBITaxa
#         from collections import defaultdict
#         import os
#         from tqdm import tqdm
#         from os.path import join as pjoin
#         from pandas import DataFrame, concat
#
#         out_file="1000_processed_reads/kaiju_table.tsv"
#
#         pat = "/crex/proj/uppstore2018116/moritz6/1000_processed_reads/"
#         samples = [s for s in os.listdir(pat) if "." not in s]
#         out_dict = {}
#
#         for s in tqdm(samples):
#             if not out_dict.get(s):
#                 out_dict[s]=defaultdict(int)
#                 with open(pjoin(pat, s, "reads", "kaiju", s+"_kaiju.out")) as handle:
#                     for l in tqdm(handle):
#                         out_dict[s][l.split()[2]]+=1
#
#         data = DataFrame.from_dict(out_dict)
#         data = data.fillna(0)
#
#         taxDb = NCBITaxa()
#         line= {i : taxDb.get_lineage(i) for i in tqdm(list(data.index))}
#         out = {i : taxDb.get_taxid_translator(l) if l else None for i,l in tqdm(line.items())}
#         tt = sum([list(l.keys()) for l in tqdm(out.values()) if l], [])
#         tt = list(set(tt))
#         tt = taxDb.get_rank(tt)
#
#         out_l = {k : {tt[k] : v for k,v in t.items()} if t else None for k,t in tqdm(out.items())}
#         taxos = DataFrame.from_dict(out_l).transpose()
#         taxos = taxos.fillna("NA")
#         all = concat([data, taxos.loc[data.index]], axis=1)
# 	all.to_csv(out_file)
#
#
# rule mash_matrix:
#     input : "1000_processed_reads/done"
#     output : "1000_processed_reads/mash_distances.txt"
#     shell :
#         """
#         dir=`dirname {input}`
#         find $dir/*/reads/mash -name "*.msh" > mashes.txt
#
#         for f in `cat mashes.txt`
#         do
#             mash dist -p 20 -l $f mashes.txt
#         done > {output}
#         """
#
# rule mash_report:
#     input : dist_file = "1000_processed_reads/mash_distances.txt",
#             metadata = "9999_metadata/metadata.csv",
#             sample_map = "9999_metadata/sample_map.csv"
#     output : out_pdf = "1000_processed_reads/mash_distances.pdf"
#     script : config['general']['scripts']['mash_report']
#
#
#
#
#
# rule library:
#     input :           "{path}/{sample}/reads/mash/{sample}.msh",
#             "{path}/{sample}/megahit/mapping/mapping_rates.txt",
#             "{path}/{sample}/reads/kaiju/{sample}_kaiju.out.summary"
#     output : "{path}/{sample}/done"
#     shell:
#         "touch {wildcards.path}/{wildcards.sample}/done"
#
#
#
# rule all_libraries :
#     input : all_dones ,
#     output : "1000_processed_reads/done"
#     shell : "touch 1000_processed_reads/done"


rule bbnorm:
    input : read1 = "{path}/fwd.fastq.gz",
             read2 = "{path}/rev.fastq.gz",
    output : read1 = "{path}/fwd.bbnorm.fastq.gz",
             read2 = "{path}/rev.bbnorm.fastq.gz",
             unp = "{path}/unp.bbnorm.fastq.gz",
             hist = "{path}/bbnorm.hist"
    params : minprob = config["read_processing"]["bbnorm"]["minprob"],
             target = config["read_processing"]["bbnorm"]["target"]
    threads : 20
    shell : """
bbnorm.sh in={input.read1} in2={input.read2} out={output.read1} out2={output.read2} threads={threads}  hist={output.hist} prefilter=t minprob={params.minprob} target={params.target}

touch {output.unp}
    """
