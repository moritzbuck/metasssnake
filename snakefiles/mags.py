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


def process_hmm_file(f) :
    domtblout_head = ["target_name", "target_accession", "target_len", "query_name" , "query_accession", "query_len" , "E-value","score-sequence" , "bias-sequence" , "nbr", "of","C-value", "I-value", "score-domain" , "bias-domain" , "hmm_from" , "hmm_to" , "ali_from" , "ali_to" , "env_from" , "env_to" , "acc" , "description_of_target"]
    data = pandas.read_csv(f, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,22))
    data = data.loc[data['E-value'] < 10e-9]
    pfams_dict = {p : [] for p in set(data['target_name'])}
    for a,b in data.iterrows():
        pfams_dict[b['target_name']] += [b['query_accession']]
    return pfams_dict


rule hmmer_all_mags :
    input : stats = "{path}/binning/{binner}/magstats.csv",
            cov_table = "{path}/filtered_assembly/mapping/map_table.tsv",
            folder = "{path}/binning/{binner}/clean_bins",
    output : pfam_sets = "{path}/binning/{binner}/pfam_sets.json",
             normed_mat = "{path}/binning/{binner}/normed_pfam_covs.csv"
    threads : 20
    run :
        import json
        from pathos.multiprocessing import ProcessingPool as Pool

        stats_sheet = pandas.read_csv(input.stats)
        bins =  list(stats_sheet['Unnamed: 0'])
        hmm_string = "hmmsearch --noali --cpu {cpus} --domtblout {out} {db} {seqs} > /dev/null"
        temp_pfam = pjoin(config['general']['temp_dir'], os.path.basename(config['mags']['pfams_db']))
        shutil.copy(config['mags']['pfams_db'], temp_pfam)

        commands = [hmm_string.format(cpus = 1, out = pjoin(input.folder, b, b + ".pfam.hmmsearch"), db = temp_pfam , seqs = pjoin(input.folder, b, b + ".faa")) for b in bins]
        pool = Pool(threads)

        def f(c):
            call(c, shell = True)
            print("runned : " + c)
            return "1"

        pool.map(f, commands)

        big_dict = {b : process_hmm_file(pjoin(input.folder, b, b + ".pfam.hmmsearch")) for b in bins}

        print("dumping PFAM dict")


        with open(output.pfam_sets, "w") :
            json.dump(big_dict)

        print("parsing gffs to link contigs to CDSs")

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


        pfam_cov_dict = { k : covs.loc[v].sum() for k,v in tqdm(all_pfams_as_ctgs.items())}
        pfam_pandas = pandas.DataFrame.from_dict(pfam_cov_dict, orient="index")
        pfam_pandas.index = [p.split(".")[0] for p in pfam_pandas.index]
        pfam_pandas.columns = [c.replace(".bam", "") for c in pfam_pandas.columns]

        with open("/home/moritz/dbs/pfam/pfam_dict.csv") as handle:
            pfam2id = {l.split()[0] : l.split()[1].split(".")[0]  for l in handle}
            id2pfam = {v : k for k,v in pfam2id.items()}

        with open("/home/moritz/dbs/pfam/sc_pfams.txt" ) as handle:
            sc_pfams = [s[:-1] for s in handle]

        norm_factor = pfam_pandas.loc[sc_pfams].mean().to_dict()
        normed_pfam_pandas = pandas.DataFrame.from_dict({ k : {kk : vv/norm_factor[k] for kk, vv in v.items()} for k, v in pfam_pandas.to_dict().items()})
        normed_pfam_pandas.to_csv(output.normed_mat)




rule annotate_all_mags :
    input : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins",
            unbinned = "{path}/{set}/assemblies/{assembler}/binning/{binner}/bins/bin-unbinned.fasta"
    output : folder = "{path}/{set}/assemblies/{assembler}/binning/{binner}/clean_bins",
             stats = "{path}/{set}/assemblies/{assembler}/binning/{binner}/magstats.csv"
    threads : 20
    run :
        bins = [f for f in os.listdir(input.folder) if f.endswith(".fasta")]
        if not os.path.exists(output.folder):
            os.makedirs(output.folder)

        prokka_line = "prokka --outdir {temp_out}/{prefix}  {meta} --force --prefix {prefix} --locustag {prefix} --cpus {threads} {bins}"
        checkm_line = "checkm lineage_wf -t {threads} -x fna {temp_out}/{prefix} {temp_out}/{prefix}/data > {temp_out}/{prefix}/checkm.txt"
        out_dict = {}

        for b in bins:
            meta = "--metagenome" if b == "bin-unbinned.fasta" else ""
            b_name = b[:-6].replace("_", "-" )
            prefix = "{set}_{assembler}_{binner}_{bin}".format(**wildcards, bin = b_name)
            prok = prokka_line.format(temp_out = config['general']['temp_dir'], meta = meta, prefix = prefix, threads = threads, bins = pjoin(input.folder, b) )

            call(prok, shell = True)
            if meta == "":
                call(checkm_line.format(threads= threads, temp_out = config['general']['temp_dir'], prefix = prefix), shell = True)
                shutil.rmtree(pjoin(config['general']['temp_dir'], prefix, "data"))

            shutil.move("{temp_out}/{prefix}".format(temp_out = config['general']['temp_dir'], prefix = prefix), output.folder)
            out_dict[prefix] = mag_stat(pjoin(output.folder, prefix), prefix)
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
