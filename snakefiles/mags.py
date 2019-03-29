from os.path import join as pjoin
import os
import uuid
from subprocess import Popen, PIPE, call
from Bio import SeqIO
from tqdm import tqdm
from pandas import DataFrame
import pandas
import hashlib

shell.prefix("module load bioinfo-tools bbmap samtools  BioPerl prokka    perl_modules; ")


def mag_stat(folder, bin_head) :
    bin_checkm = pjoin(folder , "checkm.txt")
    bin_genome= pjoin(folder ,bin_head + ".fna")
    bin_proteom = pjoin(folder ,bin_head + ".faa")
    out_dict = {}

    if "unbinned" not in bin_head:
        checkm_fields = ['Bin Id', 'Marker lineage', 'UID', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        with open(bin_checkm) as handle:
            dat = handle.readlines()
        dat = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l[:-1].split() ) } for l in dat[3:-1]}
        dat = dat[bin_head]

        out_dict['completeness'] = float(dat['Completeness'])
        out_dict['contamination'] = float(dat['Contamination'])
        out_dict['taxo:checkm'] = dat['Marker lineage']
        out_dict['strain_heterogeneity'] = float(dat['Strain heterogeneity'])
    else :
        out_dict['completeness'] = None
        out_dict['contamination'] = None
        out_dict['taxo:checkm'] = 'Unbinned'
        out_dict['strain_heterogeneity'] = None


    with open( bin_genome) as handle:
        fna = [s for s in SeqIO.parse(handle, "fasta")]
    with open( bin_proteom) as handle:
        faa = [s for s in SeqIO.parse(handle, "fasta")]

    out_dict['length'] = sum([len(s.seq) for s in fna])
    out_dict['nb_contigs'] = len(fna)
    out_dict['nb_proteins'] = len(faa)
    out_dict['coding_density'] = (3.0*sum([len(s.seq) for s in faa]))/sum([len(s.seq) for s in fna])
    out_dict['GC'] = float(sum([str(s.seq).count("G")+str(s.seq).count("C") for s in fna]))/out_dict['length']
    #out_dict['taxo:phylophlan'] = tax[bin_id]
    return out_dict

rule phylophlan :
    input : magstats = "{path}/magstats.csv",
            taxfile = "{path}/gtdbtk/gtdb.gtdbtk.tax"
    output : tree = "{path}/phylophlan/phylophlan.tree",
             tree_w_tax = "{path}/phylophlan/phylophlan.gtdbtk.tree"
    conda : "gtdbtk"
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

rule sourmash:
    input : gtdb_class = "{path}/gtdbtk/gtdb.gtdbtk.tax"
            def_db = ""
            
    run :
        mkdir sourmash
        cp gtdbtk/gtdb.gtdbtk.tax sourmash/db.tax
        grep -h -v ^identifier ~/dbs/sourmash/gtdb/gtdb_taxonomy.mod.tsv ~/dbs/sourmash/fungiii/fungiii.csv >> sourmash/db.tax
        cp ~/dbs/sourmash/all_sigs/* sourmash/sigs/
        mv clean_bins/*.sig  sourmash/sigs/

rule gtdbtk:
    input : magstats = "{path}/magstats.csv"
    output : taxfile = "{path}/gtdbtk/gtdb.gtdbtk.tax"
    conda : "gtdbtk"
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
    

def process_hmm_file(f) :
#    print("processing", f)
    domtblout_head = ["target_name", "target_accession", "target_len", "query_name" , "query_accession", "query_len" , "E-value","score-sequence" , "bias-sequence" , "nbr", "of","C-value", "I-value", "score-domain" , "bias-domain" , "hmm_from" , "hmm_to" , "ali_from" , "ali_to" , "env_from" , "env_to" , "acc" , "description_of_target"]
    data = pandas.read_csv(f, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,22), engine='python')
    data = data.loc[data['E-value'] < 10e-9]
    pfams_dict = {p : [] for p in set(data['target_name'])}
    for a,b in data.iterrows():
        pfams_dict[b['target_name']] += [b['query_accession']]
    return pfams_dict

def cov_table(wildcards):
    if "1500_" in wildcards.path or "1000_" in wildcards.path:
        return pjoin(wildcards.path,"../../filtered_assembly/mapping/map_table.tsv")
    else :
        return pjoin(wildcards.path, "FullAssembly/mapping/map_table.tsv")

rule tax_table:
    input : cov_table = cov_table,
            tax_table = "full_taxonomy.tax"
    output : normed_mat = "{path}/normed_reads_mat.csv",
             raw_mat = "{path}/reads_mats.csv"
    threads : 1
    run :
        import pandas
        cov_table = input.cov_table
        tax_table = input.tax_table
        covs = pandas.read_csv(cov_table, sep="\t", index_col = 0)
        tax_dict = pandas.read_csv(tax_table, index_col=0).to_dict()
        
        lens = covs['contigLen']
        covs = covs[[c for c in covs.columns if c.endswith(".bam")]]
        reads = covs.apply(lambda l : l*lens)
        
        print("normalize table")

        normed_reads = reads.apply(lambda l : l/sum(l), axis='index')
        normed_pfam_pandas.to_csv(output.normed_mat)



rule hmmer_table :
    input : stats = "{path}/magstats.csv",
            folder = "{path}/clean_bins",
            cov_table = cov_table,
            pfam_sets = "{path}/pfam_sets.json"
    output : normed_mat = "{path}/normed_pfam_covs.csv",
             raw_mat = "{path}/pfam_covs.csv"
    threads : 1
    run :
        with open(input.pfam_sets)  as handle:
            big_dict = json.load(handle)

        print("parsing gffs to link contigs to CDSs")

        stats_sheet = pandas.read_csv(input.stats)
        bins =  list(stats_sheet['Unnamed: 0'])

        prot2contig = {}
        for b in tqdm(bins):
            f = pjoin(input.folder, b, b + ".gff")
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
        all_pfams = set(sum(ctgs2pfam.values(),[]))
        all_pfams_as_ctgs = {p : [] for p in all_pfams}
        for k, v in tqdm(ctgs2pfam.items()):
            for vv in v:
                all_pfams_as_ctgs[vv] += [k]

        print("open coverage file")
        covs = pandas.read_csv(input.cov_table, sep="\t")
        clsts_file = pjoin(os.path.dirname(input.stats), "clusters.txt")

        if os.path.exists(clsts_file):
            print("fixing name issues between pretty bin names and names of the assembly")
            ass = os.path.dirname(input.cov_table).replace("/mapping",".fna")
            m = hashlib.md5()
            check_dict = {hashlib.md5(str(s.seq).encode('utf-8')).hexdigest() : s.id for s in tqdm(SeqIO.parse(ass, "fasta"))}

            trans_map = { s.id : check_dict[hashlib.md5(str(s.seq).encode('utf-8')).hexdigest()] for b in tqdm(bins) for s in SeqIO.parse(pjoin(input.folder, b, b + ".fna"), "fasta")}
            maps_trans = {v : k for k,v in trans_map.items()}

            prefs2 = set("_".join(s.split("_")[:-1]) for s in sum(all_pfams_as_ctgs.values(),[]))
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

        pfam_pandas.to_csv(output.raw_mat)

        print("open pfam metadata")
        with open("/home/moritz/dbs/pfam/pfam_dict.csv") as handle:
            pfam2id = {l.split()[0] : l.split()[1].split(".")[0]  for l in handle}
            id2pfam = {v : k for k,v in pfam2id.items()}

        with open("/home/moritz/dbs/pfam/sc_pfams.txt" ) as handle:
            sc_pfams = [s[:-1] for s in handle]

        print("normalize table")
        norm_factor = pfam_pandas.loc[sc_pfams].mean().to_dict()
        normed_pfam_pandas = pandas.DataFrame.from_dict({ k : {kk : vv/norm_factor[k] for kk, vv in v.items()} for k, v in pfam_pandas.to_dict().items()})
        normed_pfam_pandas.to_csv(output.normed_mat)


rule hmmer_all_mags :
    input : stats = "{path}/magstats.csv",
            folder = "{path}/clean_bins",
    output : pfam_sets = "{path}/pfam_sets.json"
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

        big_dict = {b : process_hmm_file(pjoin(input.folder, b, b + ".pfam.hmmsearch")) for b in tqdm(bins)}

        print("dumping PFAM dict")

        with open(output.pfam_sets, "w") as handle:
            json.dump(big_dict, handle)


rule annotate_all_mags :
    input : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins",
            unbinned = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins/bin-unbinned.fasta"
    output : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/clean_bins",
             stats = "{path}/{set}/assemblies/{assembler}/binning/{binner}/magstats.csv"
    threads : 16
    run :
        from pathos.multiprocessing import ProcessingPool as Pool
        bins = [f for f in os.listdir(input.folder) if f.endswith(".fasta")]
        if not os.path.exists(output.folder):
            os.makedirs(output.folder)

        prokka_line = "prokka --outdir {temp_out}/{prefix}  {meta} --force --prefix {prefix} --locustag {prefix} --cpus {threads} {bins} 2> /dev/null > /dev/null"
        checkm_line = "checkm lineage_wf -t {threads} -x fna {temp_out}/{prefix} {temp_out}/{prefix}/data > {temp_out}/{prefix}/checkm.txt 2> /dev/null"
        out_dict = {}

        pool = Pool(threads)

        def f(b):
            meta = "--metagenome" if b == "bin-unbinned.fasta" else ""
            b_name = b[:-6].replace("_", "-" )
            prefix = "{set}_{assembler}_{binner}_{bin}".format(**wildcards, bin = b_name)
            prok = prokka_line.format(temp_out = config['general']['temp_dir'], meta = meta, prefix = prefix, threads = 1, bins = pjoin(input.folder, b) )
            call(prok, shell = True)
            if meta == "":
                call(checkm_line.format(threads= 1, temp_out = config['general']['temp_dir'], prefix = prefix), shell = True)
                shutil.rmtree(pjoin(config['general']['temp_dir'], prefix, "data"))

            if os.path.exists(pjoin(output.folder,prefix)):
                shutil.rmtree(pjoin(output.folder,prefix))

            shutil.move("{temp_out}/{prefix}".format(temp_out = config['general']['temp_dir'], prefix = prefix), output.folder)


        prefix = lambda bin : "{set}_{assembler}_{binner}_{bin}".format(**wildcards, bin = bin[:-6].replace("_", "-" ))
        run_bins = [b for b in bins if not os.path.exists("{temp_out}/{prefix}/checkm.txt".format(temp_out = output.folder, prefix = prefix(b)))]
        if os.path.exists("{temp_out}/{prefix}/{prefix}.faa".format(temp_out = output.folder, prefix = prefix("bin-unbinned.fasta"))):
            run_bins.remove("bin-unbinned.fasta")
        print(run_bins)

        pool.map(f, run_bins)


        print("computing stats")
        out_dict = { prefix(b) : mag_stat(pjoin(output.folder, prefix(b)), prefix(b)) for b in bins}

        DataFrame.from_dict(out_dict, orient = 'index').to_csv(output.stats)



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



rule MAG_stats:
    input : bin_file = "{path}/{name}/{assembler}/MAGs/metabat_{name}_unbinned/{name}_unbinned.checkm",
    threads : 1
    run :


#        with open("phylophlan_tax.txt") as handle:
#            tax = [l.split()[:2] for l in handle.readlines()]
#            tax = {t[0] : ":".join([tt for tt in t[1].split(".") if "?" not in tt ]) for t in tax if "all_samples-" in t[0]}

        dat = dict([process_bin(p) for p in tqdm(bin_files)])
        DataFrame.from_dict(dat, orient = 'index').to_csv(out_file)
