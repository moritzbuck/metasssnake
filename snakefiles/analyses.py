import pandas
from subprocess import Popen, PIPE, call
import os
from tqdm import tqdm 

rule merge:
   params : completness = config['analyses']['merge']['completness'],
            contamination = config['analyses']['merge']['contamination']
   input: expand("/crex/proj/uppstore2018116/moritz6/1000_processed_reads/{name}/assemblies/megahit/binning/metabat/full_taxonomy.tax", name = samples)#, expand("/crex/proj/uppstore2018116/moritz6/1500_coasses/{name}/assemblies/megahit/binning/metabat/full_taxonomy.tax", name = [a for a in coasses if a!="Loclat"  and a != "MJ-time" ])
   output : magstats = "4500_assembly_analysis/magstats.csv",
            taxonomy = "4500_assembly_analysis/full_taxonomy.tax",
            good_mags = "4500_assembly_analysis/good_mags.txt",
            proteom = "4500_assembly_analysis/proteomics/all_proteoms.faa",
            genome = "4500_assembly_analysis/genomics/all_genomes.fna",
   run : 
       from numpy import logical_and
       from os.path import join as pjoin

       # cutoff based on https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/
       tax_files = input
       comp = params.completness
       cont = params.contamination
       mag_stats = [f.replace("full_taxonomy.tax","magstats.csv") for f in tax_files]
       
       genome_fold = pjoin(os.path.dirname(output.magstats), "genomics")
       proteom_fold = pjoin(os.path.dirname(output.magstats), "proteomics")
       
       print("merge magstats")
       levels = ['superkingdom','phylum','class','order','family','genus', "species", "strain"]
       all_stats = pandas.concat([pandas.read_csv(f) for f in tqdm(mag_stats)], sort=True)
       all_stats.index = all_stats['Unnamed: 0']
       all_tax = pandas.concat([pandas.read_csv(f, index_col=0) for f in tqdm(tax_files)], sort=True)
       all_tax = all_tax[levels]
       good_bins = all_stats.index[logical_and(all_stats.completeness > comp, all_stats.contamination*(1-(all_stats.strain_heterogeneity/100)) < cont)].to_series()
       

       all_stats.to_csv(output.magstats)
       all_tax.sort_values(by = levels).to_csv(output.taxonomy)
       good_bins.to_csv(output.good_mags)
       os.makedirs(pjoin(genome_fold, "genomes"), exist_ok = True)
       os.makedirs(pjoin(proteoms_fold, "proteoms"), exist_ok = True)

       for f in tqdm(mag_stats):
           fold = os.path.dirname(f)
           mag_list = os.listdir(pjoin(fold, "clean_bins"))
           for b in tqdm(mag_list):
               os.symlink(pjoin(fold, "clean_bins", b, b + ".fna"), pjoin(genome_fold, "genomes", b + ".fna"))
               os.symlink(pjoin(fold, "clean_bins", b, b + ".faa"), pjoin(proteoms_fold, "proteoms", b + ".faa"))
       """for f in `ls proteoms/`; do echo $f ; cat proteoms/$f >> all_proteoms.faa ; done"""

rule pre_cluster_proteom:
    params : id = config['analyses']['pre_cluster_proteom']['id'],
             length_cut = config['analyses']['pre_cluster_proteom']['length_cut'],
             temp_folder = config['general']['temp_dir']
    input : proteom = "4500_assembly_analysis/proteomics/all_proteoms.faa"
    output : clusters = "4500_assembly_analysis/proteomics/precluster/representative_seq.faa",
             groups = "4500_assembly_analysis/proteomics/precluster/clusters.tsv"
    threads : 20
    run : 
        temp_fasta = pjoin(params.temp_folder, "cdhit_temp.faa")
        print("running cd-hit")
        call("cd-hit -i {input} -o {output} -c {id} -M 0 -T {threads} -d 0 -s {cut}".format(input = input.proteom, output = temp_fasta, id = params.id, cut=params.length_cut, threads= threads), shell=True)
        clstrs = []
        record = []
        print("parsing clusters")

        with open(temp_fasta + ".clstr") as handle:
            for l in tqdm(handle):
                if l.startswith(">") :
                    clstrs += [record]
                    record = []
                else : 
                    dat = l.split()[2]
                    record += [dat[1:][:-3]]
        clstrs = clstrs[1:]
        
        print("outputing clusters")
        with open(output.groups, "w") as outp:
            outp.writelines(["Cluster_" + str(i+1) + "\t" + "\t".join(rec) + "\n" for i,rec in enumerate(clstrs)])
        print("outputing fasta")
        with open(output.clusters, "w") as outp :
            counter = 0
            with open(temp_fasta) as handle:
                for l in tqdm(handle):
                    if l.startswith(">") :
                        if counter > 0 : 
                            outp.writelines(record)
                        counter += 1
                        record = [ ">Cluster_" + str(counter) + "\n"]
                    else : 
                        record += [l]

                    
rule self_diamond :
    input : clusters = "4500_assembly_analysis/proteomics/precluster/representative_seq.faa",
    output : hits =  "4500_assembly_analysis/proteomics/precluster/selfhits.diamond"
    threads : 20 
    shell :"""
diamond makedb --db {input.clusters} --in {input.clusters}
diamond blastp --more-sensitive -p {threads} -f 6 -q {input.clusters} --db {input.clusters} -o {output.hits}
"""
