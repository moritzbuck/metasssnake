from os.path import join as pjoin
import os
from os.path import realpath
import uuid
from subprocess import Popen, PIPE, call
from Bio import SeqIO
from tqdm import tqdm
from pandas import DataFrame
import pandas
import hashlib

#shell.prefix("module load bioinfo-tools bbmap samtools  BioPerl prokka    perl_modules; ")


def process_checkm_file(file) : 
    checkm_fields = ['Bin Id', 'Marker lineage', 'UID', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+'\
, 'Completeness', 'Contamination', 'Strain heterogeneity']
    with open(file) as handle:
        dat = handle.readlines()
    dat = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l[:-1].split() ) } for l in dat[3:-1]}
    out_dict = {}
    for k, v in dat.items():
        out_dict[k] = {}
        out_dict[k]['completeness'] = float(dat[k]['Completeness'])
        out_dict[k]['contamination'] = float(dat[k]['Contamination'])
        out_dict[k]['taxo:checkm'] = dat[k]['Marker lineage']
        out_dict[k]['strain_heterogeneity'] = float(dat[k]['Strain heterogeneity'])
    return out_dict

def mag_stat(folder, bin_head) :
    bin_genome= pjoin(folder ,bin_head + ".fna")
    bin_proteom = pjoin(folder ,bin_head + ".faa")
    out_dict = {}

    with open( bin_genome) as handle:
        fna = [s for s in SeqIO.parse(handle, "fasta")]
    with open( bin_proteom) as handle:
        faa = [s for s in SeqIO.parse(handle, "fasta")]

    out_dict['length'] = sum([len(s.seq) for s in fna])
    out_dict['nb_contigs'] = len(fna)
    out_dict['nb_proteins'] = len(faa)
    out_dict['coding_density'] = (3.0*sum([len(s.seq) for s in faa]))/sum([len(s.seq) for s in fna])
    out_dict['GC'] = float(sum([str(s.seq).count("G")+str(s.seq).count("C") for s in fna]))/out_dict['length']
    return out_dict

rule phylophlan :
    input : magstats = "{path}/magstats.csv",
            taxfile = "{path}/gtdbtk/gtdb.gtdbtk.tax"
    output : tree = "{path}/phylophlan/phylophlan.tree",
             tree_w_tax = "{path}/phylophlan/phylophlan.gtdbtk.tree"
    threads : 20
    params : phylophlan_path = "/home/moritz/repos/github/phylophlan",
             phylophlan_exe = "phylophlan.py",
    run :
        temp_name = str(uuid.uuid4()).split("-")[-1]
        mag_fold = pjoin(os.path.dirname(input.magstats), "clean_bins")
        folder = pjoin(params.phylophlan_path, "input", temp_name)
        os.makedirs(folder, exist_ok=True)
        mag_dat = pandas.read_csv(input.magstats)
        cwd =  os.getcwd()

        for l in mag_dat['Unnamed: 0'][mag_dat.completeness > 10]:
            if "unbinned" not in l:
                shutil.copy(pjoin(mag_fold, l, l + ".faa"), folder)

        os.chdir(params.phylophlan_path)
        phlan_line = "python2.7 phylophlan.py --nproc {threads}  -u  {name}"
        call(phlan_line.format(threads = threads, name = temp_name), shell = True)
        shutil.copy("output/{name}/{name}.tree.nwk".format(name = temp_name), "{cwd}/{output}".format(cwd = cwd, output = output.tree))
        shutil.rmtree(pjoin("input", temp_name))
        shutil.rmtree(pjoin("output", temp_name))
        shutil.rmtree(pjoin("data", temp_name))
        os.chdir(cwd)

        with open("gtdbtk/gtdb.gtdbtk.tax") as handle:
            lines = handle.readlines()

        name_dict = {l.split()[0] : l.split()[0] + "." + l.split()[1].replace(";",".") for l in lines}

rule mash_lca:
    params : "/home/moritz/dbs/sourmash/gtdb/gtdb.lca.json"
    input : gtdb_class = "{path}/gtdbtk/gtdb.gtdbtk.tax",
            sig_list = "{path}/sourmash/sigs.list",
    output :sourmash_combo_class = "{path}/sourmash/combo.sourmash.tax",
            sourmash_local_db = "{path}/sourmash/local.lca.json"
    conda : "smash"
    shell :
        """
        folder=`dirname {output.sourmash_combo_class}`
        echo "extra db :" {params}
        cp {input.gtdb_class} $folder/db.tax
        sed -i 's/[dpcofgs]__//g' $folder/db.tax
        sourmash lca index  $folder/db.tax {output.sourmash_local_db}  $folder/good_sigs/*.sig
        sourmash lca classify --db {output.sourmash_local_db} {params} --query  $folder/sigs/*.sig > {output.sourmash_combo_class}
        """

rule merg_classes : 
    input : sourmash_combo_class = "{path}/sourmash/combo.sourmash.tax",
            gtdb_class = "{path}/gtdbtk/gtdb.gtdbtk.tax"
    output : full_class = "{path}/full_taxonomy.tax"
    threads : 1
    run : 
            import re
            sourmash_combo_class = input.sourmash_combo_class
            gtdb_class = input.gtdb_class
            full_class = output.full_class
            with open(gtdb_class) as handle : 
                gtdbtk_tax = {l.split(",")[0] : re.sub("[dpcofgs]__", '', l[:-1]) for l in handle.readlines() if not l.startswith('identifiers') }
          
            head_line = "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
  
            with open(sourmash_combo_class) as handle : 
                sourmash_tax = {l.split(",")[0] : re.sub("[dpcofgs]__", '', l[:-1]) for l in handle.readlines() if not l.startswith('ID') and "nomatch" not in l}
            confused = {k :  v.replace(",disagree", "") for k, v in sourmash_tax.items() if "disagree" in v}
            sourmash_tax = {k :  v.replace(",found", "") for k, v in sourmash_tax.items() if "found" in v}
            sourmash_tax.update(gtdbtk_tax)
            with open(full_class, "w" ) as handle:
                handle.writelines([head_line] + [s + "\n" for s in sourmash_tax.values()])


rule mash_mags:
    input : magstats = "{path}/magstats.csv"
    output : sig_list = "{path}/sourmash/sigs.list",
             mag_list = "{path}/sourmash/mags.list",
             out_folder = "{path}/sourmash/sigs/",
             good_sigs = "{path}/sourmash/good_sigs/"
#    conda : "smash"
    threads : 20
    run :
        import shutil
        from numpy import logical_and
        magstats = input.magstats
        sig_list = output.sig_list
        mag_list = output.mag_list
        folder = wildcards.path 
        out_folder = output.out_folder
        good_sigs = output.good_sigs
        os.makedirs(out_folder, exist_ok = True)
        mag_dat = pandas.read_csv(magstats)
        files = []
        for l in mag_dat['Unnamed: 0']:
            if "unbinned" not in l:
                files += [pjoin(folder, "clean_bins", l, l + ".fna")]
        with open(mag_list, "w") as handle: 
            handle.writelines([l + "\n" for l in files])
        call("sourmash compute -k 31 --scaled 10000 -p 20 `cat {maglist} `".format(maglist = mag_list), shell = True)
        
        print("moving the sigs")
        for f in tqdm(files) :
            shutil.move(os.path.basename(f) + ".sig", out_folder)
        
        files = [pjoin(out_folder, os.path.basename(l)) + ".sig" for l in files]
        
        def correct(f):                                                                  
            with open(f) as handle:                                                                                                                          
                tt = json.load(handle)                                                                                                            
                tt[0]['name'] = tt[0]['filename'].split("/")[-1].replace(".fna","")
            with open(f,"w") as handle:
                json.dump(tt,handle)
        
        print("correcting the names of the sigs")
        for f in tqdm(files):
              correct(f)

        print("moving sigs of good sigs")
        goods =  mag_dat.loc[logical_and(mag_dat.completeness > 30 , mag_dat.contamination < 5)]
        for l in goods['Unnamed: 0']:
            if "unbinned" not in l:
                shutil.move(pjoin(out_folder, os.path.basename(l) + ".fna.sig"), good_sigs )
        

        with open(sig_list, "w") as handle: 
            handle.writelines([l + "\n" for l in files])


rule gtdbtk:
    input : magstats = "{path}/magstats.csv"
    output : taxfile = "{path}/gtdbtk/gtdb.gtdbtk.tax"
    threads : 20
    run :
        folder = os.path.dirname(output.taxfile)
        os.makedirs(folder, exist_ok=True)
        mag_dat = pandas.read_csv(input.magstats)
        for l in mag_dat['Unnamed: 0']:
            if "unbinned" not in l:
                shutil.copy(pjoin(folder, "../clean_bins", l, l + ".fna"), folder)

        identify_line = "gtdbtk identify --genome_dir {path} --out_dir {path}/temp -x fna --cpus {threads}"
        align_line = "gtdbtk align --identify_dir {path}/temp --out_dir {path}/temp --cpus {threads} "
        classify_line = "gtdbtk classify --genome_dir {path} --align_dir {path}/temp --out_dir {path}/temp --cpus 1"

        call(identify_line.format(threads= threads, path = folder), shell = True)
        call(align_line.format(threads= threads, path = folder), shell = True)
        call(classify_line.format(threads= threads, path = folder), shell = True)

        for f in os.listdir(folder):
            if f.endswith(".fna"):
                os.remove(pjoin(folder, f))

        head_line = "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
        with open(pjoin(folder, "temp", "gtdbtk.ar122.summary.tsv")) as handle:
            arc_data = [l.split('\t')[0:2] for l in handle.readlines()[1:]]
            arc_data = [t[0] + "," + t[1].replace(";",",") + ",\n"  for t in arc_data]
        with open(pjoin(folder, "temp", "gtdbtk.bac120.summary.tsv")) as handle:
            bac_data = [l.split('\t')[0:2] for l in handle.readlines()[1:]]
            bac_data = [t[0] + "," + t[1].replace(";",",") + ",\n"  for t in bac_data]
        with open(output.taxfile, "w") as handle:
            handle.writelines([head_line] + arc_data + bac_data)
    

def process_hmm_file(f, cutoff) :
#    print("processing", f)
    domtblout_head = ["target_name", "target_accession", "target_len", "query_name" , "query_accession", "query_len" , "E-value","score-sequence" , "bias-sequence" , "nbr", "of","C-value", "I-value", "score-domain" , "bias-domain" , "hmm_from" , "hmm_to" , "ali_from" , "ali_to" , "env_from" , "env_to" , "acc" , "description_of_target"]
    with open(f) as handle : 
        dd = {i : { k: ff for k, ff in zip(domtblout_head, f.split()[:len(domtblout_head)])} for i,f in enumerate(handle) if not f.startswith("#")}
    if len(dd) == 0 :
        return ({},pandas.DataFrame())
    data = pandas.DataFrame.from_dict(dd, orient="index")
#    data = pandas.read_csv(f, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,22), engine='python')
#    if len(data) > 0 and data['target_accession'][0] != '-':
#        data = pandas.read_csv(f, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,22))
    pfams_dict = {p : [] for p in set(data['target_name'])}
    for a,b in data.loc[data['I-value'].apply(float) < cutoff].iterrows():
        pfams_dict[b['target_name']] += [b['query_accession']]
    return (pfams_dict, data)

def cov_table(wildcards):
    if "1500_" in wildcards.path or "1000_" in wildcards.path:
        return realpath(pjoin(wildcards.path,"../../filtered_assembly/mapping/map_table.tsv"))
    else :
        return pjoin(wildcards.path, "FullAssembly/mapping/map_table.tsv")

rule tax_table:
    input : cov_table = cov_table,
            tax_table = "{path}/full_taxonomy.tax",
            stats = "{path}/magstats.csv",
    output : bin_wise = "{path}/abundance_tables/abundance_per_bin.csv",
             folder = "{path}/abundance_tables/"
    threads : 1
    run :
        import pandas
        from pandas import Series
        import os 

        cov_table = input.cov_table
        tax_table = input.tax_table
        stats = input.stats
        folder = pjoin(wildcards.path, "clean_bins")
        bin_wise = output.bin_wise
        out_folder = output.folder

        covs = pandas.read_csv(cov_table, sep="\t", index_col = 0)
        stats_sheet = pandas.read_csv(stats)
        bins =  list(stats_sheet['Unnamed: 0'])

        clsts_file = pjoin(os.path.dirname(stats), "clusters.txt")

        levels = ['superkingdom','phylum','class','order','family','genus', "species", "strain"]
        tax_df = pandas.read_csv(tax_table, index_col=0).fillna("")
        clean_tax =  lambda s: Series([ss if ss != "" else ("undetermined" + last([tt for tt in s if tt ])) for ss in s], levels )
        last = lambda ll : ("_" + ll[-1]) if len(ll) > 0 else ""

        tax_df = tax_df.apply(clean_tax, axis=1)


        if os.path.exists(clsts_file):
            print("fixing name issues between pretty bin names and names of the assembly")
            ass = os.path.dirname(cov_table).replace("/mapping",".fna")
            m = hashlib.md5()
            check_dict = {hashlib.md5(str(s.seq).encode('utf-8')).hexdigest() : s.id for s in tqdm(SeqIO.parse(ass, "fasta"))}

            trans_map = { s.id : check_dict[hashlib.md5(str(s.seq).encode('utf-8')).hexdigest()] for b in tqdm(bins) for s in SeqIO.parse(pjoin(folder, b, b + ".fna"), "fasta")}
            maps_trans = {v : k for k,v in trans_map.items()}

            prefs2 = set("_".join(s.split("_")[:-1]) for s in covs.index)
            assert len(prefs2) == 1
            prefs2 = list(prefs2)[0] + "_bin-"

            covs.index = [maps_trans.get(c,c).replace("bin-",prefs2) for c in covs.index]
            covs = covs.loc[[c for c in covs.index if "/" in c]]
            tax_df.index = [t.split("_")[-1].replace("bin-",prefs2) for t in tax_df.index]
        else :
            covs.index = covs.index



        unbinned = Series(['unbinned']*len(levels), index =levels)
        unclassified = Series(['undetermined']*len(levels), index =levels)

        lens = covs['contigLen']
        covs = covs[[c for c in covs.columns if c.endswith(".bam")]]
        reads = covs.apply(lambda l : l*lens)
        reads.columns = [c.replace(".bam", "") for c in reads.columns]
        print("normalize table")
        def get_tax(mab) :
            if mab in tax_df.index :
                return tax_df.loc[mab]
            elif 'unbinned' in mab:
                return unbinned
            else :
                return unclassified
     
        normed_reads = reads.apply(lambda l : l/sum(l), axis='index')
        normed_reads['contig'] = normed_reads.index
        normed_reads['mab'] = [c.split(":")[0] for c in normed_reads.contig]
        normed_reads = normed_reads.apply(lambda l : l.append( get_tax(l.mab)) , axis=1)

        melted_reads = normed_reads.melt(id_vars = ['contig', 'mab'] + levels, value_name="reads", var_name="sample")
        
        value_dict = {g[0] : g[1].reads.sum() for g in tqdm(melted_reads.groupby(['mab', 'sample']))}
        contig_wise_table = pandas.DataFrame.from_dict({ sample : { k[0] : v for k, v in value_dict.items() if k[1] == sample}  for sample in  tqdm(normed_reads.columns) if sample not in ['contig'] + levels })
        contig_wise_table['mab'] = contig_wise_table.index        
        contig_wise_table = contig_wise_table.apply(lambda l : l.append( get_tax(l.mab)) , axis=1)


        melted_contig_read = contig_wise_table.melt(id_vars = ['mab'] + levels , value_name="reads", var_name="sample")

        level_tables = {}
        for l in levels:
            value_dict = {g[0] : g[1].reads.sum() for g in melted_contig_read.groupby([l, "sample"])}
            level_tables[l] = pandas.DataFrame.from_dict({ sample : { k[0] : v for k, v in value_dict.items() if k[1] == sample}  for sample in  normed_reads.columns if sample not in ['contig'] + levels })
            
        contig_wise_table.to_csv(bin_wise)
        reads.to_csv(pjoin(out_folder, "raw_read_counts.csv"))
        normed_reads.to_csv(pjoin(out_folder,"abundance_per_contig.csv"))
        for k, l in level_tables.items():
            l.to_csv(pjoin(out_folder,"abundance_per_" + k + ".csv"))




rule hmmer_table :
    input : stats = "{path}/magstats.csv",
            folder = "{path}/clean_bins",
            cov_table = cov_table,
            pfam_sets = "{path}/pfam_sets.json"
    output : normed_mat = "{path}/normed_pfam_covs.csv",
             raw_mat = "{path}/pfam_covs.csv"
    threads : 1
    run :
        import json
        folder = input.folder 
        cov_table = input.cov_table
        normed_mat = output.normed_mat
        raw_mat = output.raw_mat
        pfam_sets = input.pfam_sets
        with open(pfam_sets)  as handle:
            big_dict = json.load(handle)


        print("parsing gffs to link contigs to CDSs")

        stats_sheet = pandas.read_csv(input.stats)
        bins =  list(stats_sheet['Unnamed: 0'])
        prot2contig = {}
        for b in tqdm(bins):
            f = pjoin(folder, b, b + ".gff")
            with open(f) as handle:
                for line in handle:
                    if line == '##FASTA\n':
                        break
                    elif line[0] != "#":
                        cont_id = line.split("\t")[0]
                        traits = {t.split("=")[0] : t.split("=")[1]  for t in line.split("\t")[-1].split(";")}
                        cds_id = traits.get('ID')
                        prot2contig[cds_id] = cont_id

        ctgs2pfam = {g + ":" +  prot2contig[k].split(":")[1] : v  for g,d in tqdm(big_dict.items()) for k, v in d.items()}
        all_pfams = set()
        for g,d in tqdm(big_dict.items()):
            for k, v in d.items():
                all_pfams = all_pfams.union(v)
#        all_pfams = set(sum(ctgs2pfam.values(),[]))
        all_pfams_as_ctgs = {p : [] for p in all_pfams}
        for k, v in tqdm(ctgs2pfam.items()):
            for vv in v:
                all_pfams_as_ctgs[vv] += [k]

        print("open coverage file")
        covs = pandas.read_csv(input.cov_table, sep="\t")
        clsts_file = pjoin(os.path.dirname(input.stats), "clusters.txt")

        if os.path.exists(clsts_file):
            print("fixing name issues between pretty bin names and names of the assembly")
            ass = os.path.dirname(cov_table).replace("/mapping",".fna")
            m = hashlib.md5()
            check_dict = {hashlib.md5(str(s.seq).encode('utf-8')).hexdigest() : s.id for s in tqdm(SeqIO.parse(ass, "fasta"))}

            trans_map = { s.id : check_dict[hashlib.md5(str(s.seq).encode('utf-8')).hexdigest()] for b in tqdm(bins) for s in SeqIO.parse(pjoin(folder, b, b + ".fna"), "fasta")}
            maps_trans = {v : k for k,v in trans_map.items()}

            prefs2 = set("_".join(s.split("_")[:-1]) for s in ctgs2pfam)
            assert len(prefs2) == 1
            prefs2 = list(prefs2)[0] + "_bin-"

            covs.index = [maps_trans.get(c,c).replace("bin-",prefs2) for c in covs['contigName']]
            covs = covs.loc[[c for c in covs.index if "/" in c]]
        else :
            covs.index = covs['contigName']

        lengDict = covs['contigLen'].to_dict()
        covs = covs[[c for c in covs.columns if c.endswith(".bam")]]


        print("make coverage dict")
        pfam_cov_dict = { k : covs.loc[v].sum() for k,v in tqdm(all_pfams_as_ctgs.items())}
        pfam_pandas = pandas.DataFrame.from_dict(pfam_cov_dict, orient="index")
        pfam_pandas.index = [p.split(".")[0] for p in pfam_pandas.index]
        pfam_pandas.columns = [c.replace(".bam", "") for c in pfam_pandas.columns]

        pfam_pandas.to_csv(raw_mat)

        print("open pfam metadata")
        with open("/home/moritz/dbs/pfam/pfam_dict.csv") as handle:
            pfam2id = {l.split()[0] : l.split()[1].split(".")[0]  for l in handle}
            id2pfam = {v : k for k,v in pfam2id.items()}

        with open("/home/moritz/dbs/pfam/sc_pfams.txt" ) as handle:
            sc_pfams = [s[:-1] for s in handle]

        print("normalize table")
        norm_factor = pfam_pandas.loc[sc_pfams].mean().to_dict()
        normed_pfam_pandas = pandas.DataFrame.from_dict({ k : {kk : vv/norm_factor[k] for kk, vv in v.items()} for k, v in pfam_pandas.to_dict().items()})
        normed_pfam_pandas.to_csv(normed_mat)


rule hmmer_all_mags :
    input : stats = "{path}/magstats.csv",
            folder = "{path}/clean_bins",
    output : pfam_sets = "{path}/pfam_sets.json",
             raw_hits = "{path}/raw_pfam_hits.csv"
    threads : 20
    run :
        import json
        from pathos.multiprocessing import ProcessingPool as Pool

        stats_sheet = pandas.read_csv(input.stats)
        bins =  list(stats_sheet['Unnamed: 0'])
        hmm_string = "hmmsearch --noali --cpu {cpus} --domtblout {out} {db} {seqs} > /dev/null"
        temp_pfam = pjoin(config['general']['temp_dir'], os.path.basename(config['mags']['pfams_db']))
        shutil.copy(config['mags']['pfams_db'], temp_pfam)
        bad_bins = set([b for b in bins if os.stat(pjoin(input.folder, b, b + ".faa")).st_size == 0])
        for b in bad_bins :
            call("touch " + pjoin(input.folder, b, b + ".pfam.hmmsearch"), shell = True)
        commands = [hmm_string.format(cpus = 1, out = pjoin(input.folder, b, b + ".pfam.hmmsearch"), db = temp_pfam , seqs = pjoin(input.folder, b, b + ".faa")) for b in bins if not os.path.exists(pjoin(input.folder, b, b + ".pfam.hmmsearch"))]
        print(commands)
        pool = Pool(threads)

        def f(c):
            call(c, shell = True)
            return "1"


        pool.map(f, commands)

        big_dict = {b : process_hmm_file(pjoin(input.folder, b, b + ".pfam.hmmsearch"), config['mags']['pfam_cutoff']) for b in tqdm(bins)}
        all_raw_hits = pandas.concat([ v[1]  for v in big_dict.values()], axis = 0)
        big_dict = {k : v[0] for k, v in big_dict.items()}
        print("dumping PFAM dict")
        all_raw_hits.to_csv(output.raw_hits)
        with open(output.pfam_sets, "w") as handle:
            json.dump(big_dict, handle)

rule prokka_all: 
    input : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins",
            unbinned = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins/bin-unbinned.fasta"
    output : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/clean_bins",
            flag = "{path}/{set}/assemblies/{assembler}/binning/{binner}/.all_annotated"
    threads : 20 
    run : 
        from pathos.multiprocessing import ProcessingPool as Pool
        bins = [f for f in os.listdir(input.folder) if f.endswith(".fasta")]
        if not os.path.exists(output.folder):
            os.makedirs(output.folder)

        prokka_line = "prokka --outdir {temp_out}/{prefix}  {meta} --force --prefix {prefix} --locustag {prefix} --cpus {threads} {bins} 2> /dev/null > /dev/null; echo {prefix} done"

        pool = Pool(threads)

        def f(b):
            meta = "--metagenome" if b == "bin-unbinned.fasta" else ""
            b_name = b[:-6].replace("_", "-" )
            prefix = "{set}_{assembler}_{binner}_{bin}".format(**wildcards, bin = b_name)
            prok = prokka_line.format(temp_out = config['general']['temp_dir'], meta = meta, prefix = prefix, threads = 1, bins = pjoin(input.folder, b) )
            call(prok, shell = True)

            if os.path.exists(pjoin(output.folder,prefix)):
                shutil.rmtree(pjoin(output.folder,prefix))

            shutil.move("{temp_out}/{prefix}".format(temp_out = config['general']['temp_dir'], prefix = prefix), output.folder)


        prefix = lambda bin : "{set}_{assembler}_{binner}_{bin}".format(**wildcards, bin = bin[:-6].replace("_", "-" ))

        run_bins = [b for b in bins if not os.path.exists("{temp_out}/{prefix}/{prefix}.faa".format(temp_out = output.folder, prefix = prefix(b)))]

        if "bin-unbinned.fasta" in run_bins:
            run_bins.remove("bin-unbinned.fasta")
            run_bins = ["bin-unbinned.fasta"] + run_bins

        print(len(run_bins), " bins to predict with prokka")

        pool.map(f, run_bins)
        call("touch "  +  output.flag, shell=True)

rule checkm_all :
    input : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/clean_bins",
            flag = "{path}/{set}/assemblies/{assembler}/binning/{binner}/.all_annotated"
    output : chckm = "{path}/{set}/assemblies/{assembler}/binning/{binner}/checkm.txt"
    threads : 20 
    run :
        bins = [f for f in os.listdir(input.folder)]

        temp_folder = pjoin(config['general']['temp_dir'], "checkm_" + wildcards.set)
        print("copying bins to temp" )
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

        for b in tqdm(bins):
            if "unbinned" not in b:
                shutil.copy(pjoin(input.folder, b, b +".fna"), temp_folder)

        checkm_line = "checkm lineage_wf -t {threads} -x fna {folder} {folder}/data > {outfile} "

        call(checkm_line.format(threads= threads, folder = temp_folder, outfile = output.chckm), shell = True)       

rule MAG_stats :
    input : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/clean_bins",
            flag = "{path}/{set}/assemblies/{assembler}/binning/{binner}/.all_annotated",
            checkm = "{path}/{set}/assemblies/{assembler}/binning/{binner}/checkm.txt"
    output :stats = "{path}/{set}/assemblies/{assembler}/binning/{binner}/magstats.csv"
    threads : 1
    run :
        from tqdm import tqdm 
        bins = [f for f in os.listdir(input.folder)]

        print(len(bins), " bins to process")

        print("computing basic MAG stats")
        stat_dict = { b : mag_stat(pjoin(input.folder, b), b) for b in tqdm(bins)}

        print("parsing checkm data")
        checkm_dict = process_checkm_file(input.checkm)

        print("parsing taxonomy")
#        with open(input.full_class, "w" ) as handle:
#            handle.readlines([head_line] + [s + "\n" for s in sourmash_tax.values()])
#            head_line = handle.readline()[:-1].split()[1:]
#            taxo_dict = {l.split()[0] : {k :v for k,v in zip(head_line , l.split(1:))} for l in handle}

        for k, v in stat_dict.items():
            if checkm_dict.get(k):
                v.update(checkm_dict.get(k))
#            if taxo_dict.get(k):
#            v.update(taxo_dict.get(k))

        DataFrame.from_dict(stat_dict, orient = 'index').to_csv(output.stats)



rule filter_good_MAGs :
# cutoff based on https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/
    params : contamination = 5, completeness = 40
    input : path = "{path}/{name}/bins/{assembler}/MAGs",
            anots = "{path}/{name}/bins/{assembler}/annotated",
            checkm = "{path}/{name}/bins/{assembler}/MAGs/{name}.checkm"
    output : "{path}/{name}/bins/{assembler}/good_MAGs"
    run :
        from os.path import join as pjoin
        import os
        import shutil

        checkm_fields = ['Bin Id', 'Marker lineage', 'UID', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        with open(input.checkm) as handle:
            dat = handle.readlines()
        dat = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l[:-1].split() ) } for l in dat[3:-1]}
        good_bact_MAGs = [k for k,v in dat.items() if float(v['Completeness']) > params.completeness and float(v['Contamination']) < params.contamination]

        os.makedirs(output[0])
        for f in good_bact_MAGs:
            shutil.copytree(pjoin(input.path, f), pjoin(output[0],f))



