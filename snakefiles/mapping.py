from subprocess import call
from os.path import join as pjoin
import os
import shutil
from subprocess import Popen, PIPE

def all_samples(wildcards):
        import pandas
        if config['general'].get('exp_json'):
            with open(config['general'].get('exp_json')) as handle:
                sets = json.load(handle)
                sample = [ f for
f in wildcards.path.split("/") if f in  sum(sets.values(),[]) ]
                assert len(sample) < 2
                if len(sample) == 1:
                    sety = [k for k, v in sets.items() if sample[0] in v]
                    assert len(sety)==1
                    sample_from_sets = sets[sety[0]]
            return sample_from_sets
        else :
            path = "1000_processed_reads/"
            samples = [d for d in os.listdir(path) if os.path.isdir(pjoin(path,d)) ]
            return samples

def all_bams(wildcards):
    samples = all_samples(wildcards)
    path = "{path}/mapping/bams/".format(path = wildcards.path)
    return [pjoin(path,s + ".bam") for s in samples]


def all_clean_libs(wildcards):
    files = os.listdir(config['general']['raw_folder'])
    fwds = list({"1000_processed_reads/{sample}/reads/fwd.fastq.gz".format(sample = f.split("_")[0]) for f in files })
    return fwds

rule bbmap_index:
    input : ref = "{path}/{fasta}.fna"
    output : folder_ref = "{path}/{fasta}/mapping/ref",
             folder = "{path}/{fasta}/mapping",
             gz = "{path}/{fasta}/mapping/ref/genome/1/chr1.chrom.gz",
             flag = "{path}/{fasta}/mapping/ref.ed"
    run :
       call("bbmap.sh ref={ref} path={folder}".format(ref = input.ref, folder = output.folder), shell = True)
       call("touch " + output.flag, shell = True)

rule sample_wise_bbmap :
    input : index = "{path}/mapping/ref.ed",
            ref_path = "{path}/mapping",
            fwd = "1000_processed_reads/{sample}/reads/fwd.fastq.gz",
            rev = "1000_processed_reads/{sample}/reads/rev.fastq.gz",
    output : bam = "{path}/mapping/bams/{sample}.bam",
             wdups_stats = "{path}/mapping/bams/{sample}_sorted.stats",
             stats = "{path}/mapping/bams/{sample}.stats",
    threads : 20
    run :
        bb_string = "bbmap.sh  in={fwd} in2={rev} threads={threads} out={out} bamscript={bams} path={ref}"
        temp_bam = pjoin(config['general']['temp_dir'], wildcards.sample + ".sam")
        bamsc = pjoin(config['general']['temp_dir'], "bamscr.sh")
        call(bb_string.format(fwd = input.fwd, rev = input.rev, threads = threads, ref = input.ref_path, out = temp_bam, bams = bamsc), shell = True, stderr = PIPE)
        call(bamsc, shell=True)
        call("sambamba flagstat -t {threads} {tdir}/{samp}_sorted.bam > {wdup}".format(threads = threads, samp = wildcards.sample, wdup = output.wdups_stats, tdir = config['general']['temp_dir']), shell=True)
        call("samtools rmdup  {tdir}/{samp}_sorted.bam {tdir}/{samp}.bam 2> /dev/null".format(samp = wildcards.sample, tdir = config['general']['temp_dir']), shell = True)
        call("samtools index {tdir}/{sample}.bam". format(sample = wildcards.sample, tdir = config['general']['temp_dir']), shell = True)
        call("sambamba flagstat  -t {threads} {tdir}/{samp}.bam > {stats}".format(threads = threads, samp = wildcards.sample, stats = output.stats, tdir = config['general']['temp_dir']), shell = True)
        for f in os.listdir(config['general']['temp_dir']):
            if f.startswith(wildcards.sample + ".bam"):
                  shutil.move(pjoin(config['general']['temp_dir'], f), os.path.dirname(output.bam))


rule bbmap_all_samples:
    input : "{path}/mapping/ref.ed", all_bams
    output : "{path}/mapping/map_table.tsv", "{path}/mapping/paired_contigs.tsv"
    threads : 20
    shell : """
    jgi_summarize_bam_contig_depths --outputDepth {output[0]}  --pairedContigs {output[1]}  `dirname {input[1]}`/*.bam
    """

rule bbmap_bining_map:
    input : diag_map = "{path}/mapping/mapping_rates.txt"
    output : sample_list = dynamic("{path}/mapping/bams/{sample}.ph")
    threads : 1
    run :
        import pandas
        rates_dict = {k : v['raw_map'] for k,v in pandas.read_csv(input.diag_map, index_col=0).transpose().to_dict()
.items()}
        if config['general'].get('exp_json'):
            with open(config['base'].get('exp_json')) as handle:
                sets = json.load(handle)
                sample = [ f for
f in wildcards.path.split("/") if f in  sum(sets.values(),[]) ]
                assert len(sample) < 2
                if len(sample) == 1:
                    sety = [k for k, v in sets.items() if sample[0] in v]
                    assert len(sety)==1
                    sample_from_sets = sets[sety[0]]
            rates_dict = {k : v for k, v in rates_dict.items() if k in sample_from_sets}

        vvs = sorted(rates_dict.items(),reverse = True, key = lambda i: i[1])
        vvs = vvs[:config['binning']['bin_mapping_libs']]
        vvs = [v for v in vvs if v[1] > config['binning']['bin_map_min']]
        for v in vvs:
            call("touch {path}/mapping/bams/{sample}.ph".format(path = wildcards.path, sample=v) , shell = True)
#        with open(output.sample_list, "w") as handle:
#            handle.writelines([ v[0] + "," + v[1] + "\n" for v in vvs ])

rule bbmap_diagnostic:
    input : ref_path = "{path}/mapping",
            libs = all_clean_libs,
            refd = "{path}/mapping/ref.ed"
    output : table = "{path}/mapping/mapping_rates.txt"
    threads : 20
    run :
        import pandas
        out_dict = {}
        for fwd in input.libs:
            sample = fwd.split("/")[1]
            rev = fwd.replace("fwd.fastq.gz", "rev.fastq.gz")
            print("echo mapping {sample} to {ref}".format(sample = sample, ref = input.ref_path))
            bb_string = "bbmap.sh  in={fwd} in2={rev} threads={threads} out=/dev/null reads={reads} path={ref}"
            process = Popen(bb_string.format(fwd = fwd, rev = rev, threads = threads, ref = input.ref_path, reads = config['mapping']['diagnostic']['reads']), shell = True, stderr = PIPE)
            out, err = process.communicate()
            out_dat = [l for l in err.decode().split("\n") ]
            raw_map = sum([float(l.split()[1][:-1]) for l in out_dat if l.startswith("mapped:")])/2
            mated_map = [float(l.split()[2][:-1]) for l in out_dat if l.startswith("mated pairs: ")][0]
            out_dict[sample] = {
            'raw_map' : raw_map,
            'mated_map' : mated_map
            }
        pandas.DataFrame.from_dict(out_dict, orient = "index").to_csv(output.table)
