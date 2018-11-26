from subprocess import call
from os.path import join as pjoin
import os
import shutil

shell.prefix("module load bioinfo-tools bbmap samtools; ")

def all_samples():
    path = "1000_processed_reads/"
    samples = [d for d in os.listdir(path) if os.path.isdir(pjoin(path,d)) ]
    return samples


def all_bams(wildcards):
    samples = all_samples()
    path = "{path}/mapping/bams/".format(path = wildcards.path)
    return [pjoin(path,s + ".bam") for s in samples]


def all_bin_samples(wildcards):
    path = "{path}/mapping/".format(path = wildcards.path)
    folds = [i for i,v in  enumerate(path.split("/")) if v == "1500_assemblies"]
    if len(folds) == 1:
        coas_name = path.split("/")[folds[0]+1]
        with open(pjoin(COAS_FILES_DIR, coas_name + ".txt")) as handle:
            samples_from_coas = [l.strip() for l in handle]
    else :
        samples_from_coas = []
    rates_file = pjoin(path,"mapping_rates.txt")
    with open( rates_file ) as  handle:
        rates = {l.split()[0] : float(l.split()[1]) for l in handle.readlines()[1:] if float(l.split()[1]) == float(l.split()[1])}
    vvs = sorted(list(rates.values()),reverse = True)
    cutoff = vvs[-1] if len(vvs) < BIN_MAPPING_LIBS else vvs[BIN_MAPPING_LIBS]
    cutoff = cutoff if cutoff > BIN_MAP_MIN else BIN_MAP_MIN
    samples = [k for k, v in rates.items() if v > cutoff]
    samples = list(set(samples_from_coas).union(samples))
    if CLEAN_BINNING and len(samples_from_coas) > 1:
        samples = samples_from_coas
    return [pjoin(pjoin(path,"bams",s + ".bam")) for s in samples]

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

rule samnple_wise_bbmap :
    input : index = "{path}/mapping/ref/genome/1/chr1.chrom.gz",
            fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
            rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
    output : bam = "{path}/mapping/bams/{sample}.bam",
#             wdups_stats = "{path}/mapping/bams/{sample}_sorted.stats",
#             stats = "{path}/mapping/bams/{sample}.stats"
    threads : 20
    shell : """
        module load bioinfo-tools
        module load bbmap
        module load samtools

        home=`pwd`
        cd `dirname {input.index}`/../../../

        bbmap.sh  in=$home/{input.fwd} in2=$home/{input.rev} threads={threads} bamscript=/scratch/{wildcards.sample}.sh out=/scratch/{wildcards.sample}.sam

        /scratch/{wildcards.sample}.sh
        sambamba flagstat -t {threads} /scratch/{wildcards.sample}_sorted.bam > $home/{wildcards.sample}_sorted.stats
        samtools rmdup  /scratch/{wildcards.sample}_sorted.bam /scratch/{wildcards.sample}.bam
        samtools index /scratch/{wildcards.sample}.bam
        sambamba flagstat  -t {threads} /scratch/{wildcards.sample}.bam > $home/{wildcards.sample}.stats

        rm /scratch/{wildcards.sample}.sam
        rm /scratch/{wildcards.sample}_sorted.bam*
        mv /scratch/{wildcards.sample}.bam* bams/
    """


#rule summrize_bbmap:
#    input

rule bbmap_all_samples:
    input : "{path}/mapping/ref/genome/1/chr1.chrom.gz", all_bams
    output : "{path}/mapping/map_table.tsv",
    threads : 1
    shell : """
    jgi_summarize_bam_contig_depths --outputDepth {input[0]}/map_table.tsv  --pairedContigs {input[0]}/paired_contigs.tsv  {input[0]}/bams/*.bam
    """

rule bbmap_binning_samples:
    params : home = pjoin(config['general']['home'],"1000_processed_reads"),
    input : index = "{path}/mapping/ref/genome/1/chr1.chrom.gz",
    	    path = "{path}/mapping/",
            bams = all_bin_samples
    output : "{path}/mapping/binmap_table.tsv",
    threads : 1
    shell : """
    jgi_summarize_bam_contig_depths --outputDepth {input.path}/binmap_table.tsv  --pairedContigs {input.path}/binpaired_contigs.tsv {input.bams}
    """


rule bbmap_diagnostic:
    params : reads = config['mapping']['diagnostic']['reads']
    input : ref = "{path}/mapping/ref.ed", 
            libs = all_clean_libs
    output : "{path}/mapping/mapping_rates.txt",
    threads : 20
    run : 
       for fwd in input.libs:
           sample = fwd.split("/")[1]
           rev = fwd.replace("fwd.fastq.gz", "rev.fastq.gz")
           print("echo mapping {sample} to {ref}".format(sample = sample, ref = input.ref))
           bb_string = "bbmap.sh  in={fwd} in2={rev} threads={threads} out=/dev/null reads={reads}"
           call(bb_string.format(fwd = fwd, rev = rev, threads = threads), shell = True, 
           cat tmp | grep -E "mated|mapped" | cut -f2  | tr -d ' ' | tr -d % | tr '\n'  '\t' >>  $home/{output}
           echo >> $home/{output}
           rm tmp

       pass
       """
        module load bioinfo-tools
        module load bbmap
        module load samtools

        home=`pwd`
        cd {wildcards.path}/mapping

        echo $'sample\tmated\tfwd\trev' > $home/{output}

        for s in `echo """ + " ".join(all_samples()) + """ | tr ' ' "\n"`
        do
            echo mapping $s to {wildcards.path}
            base={params.home}/$s/reads/trimmomatic/$s
            bbmap.sh  in=${{base}}_1P.fastq.gz in2=${{base}}_2P.fastq.gz threads={threads} out=/dev/null reads={params.reads} 2> tmp
            
echo -n $s $'\t' >>  $home/{output}
            cat tmp | grep -E "mated|mapped" | cut -f2  | tr -d ' ' | tr -d % | tr '\n'  '\t' >>  $home/{output}
            echo >> $home/{output}
            rm tmp
        done
    """
