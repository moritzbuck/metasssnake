import pandas
from subprocess import Popen, PIPE, call
import os
from tqdm import tqdm

rule phigaro :
   params : prop = 0.25
   input : genome = "{path}/genomics/all_genomes.fna",
           taxonomy = "{path}/full_taxonomy.tax",
           magstats = "{path}/magstats.csv",
   output : tax_file = "{path}/full_taxonomy_w_virus.tax",
            phigaro_file = "{path}/genomics/phigaro.txt"
   threads : 20
   run : 
      metagenome = input.genome
      ori_tax = input.taxonomy
      new_tax = output.tax_file
      magstats_file = input.magstats
      outfile = output.phigaro_file
      prop = params.prop

#      call("phigaro -f {genome} -t 20 -o {ouput}".format(genome = metagenome, output = outfile) , shell = True  )
      call("cp 2000_MAG_sets/Loclat/genomics/phigaro.txt.bk 2000_MAG_sets/Loclat/genomics/phigaro.txt", shell = True)
      with open(outfile) as handle:
          phi = [(l.split()[0][1:].split(":")[0], float(l.split()[2])-float(l.split()[1])+1) for l in handle if not l.endswith("end\n")]
      mag2phicov = {p[0] : 0  for p in phi} 
      for p in phi:
          mag2phicov[p[0]] += p[1]

      magstats = pandas.read_csv(magstats_file, index_col=0)
      magstats['vir_fract'] = [mag2phicov.get(k,0) for k in magstats.index]
      magstats.to_csv(magstats_file)
      viruses = magstats.loc[magstats.vir_fract > magstats.length*prop].index
      with open(ori_tax) as handle:
          tax = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}
      
      for k in viruses : 
          tax[k] = ['Virus']
      
      del tax['identifiers']
      head_line = "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
      with open(new_tax, "w" ) as handle:
          handle.writelines([head_line] + [",".join([k] + v) + "\n" for k,v in tax.items()])



rule merge:
   params : completness = config['analyses']['merge']['completness'],
            contamination = config['analyses']['merge']['contamination']
   input: expand("1000_processed_reads/{name}/assemblies/megahit/binning/metabat/full_taxonomy.tax", name = samples), expand("1500_coasses/{name}/assemblies/megahit/binning/metabat/full_taxonomy.tax", name = coasses)
   output : magstats = "4500_assembly_analysis/magstats.csv",
            taxonomy = "4500_assembly_analysis/full_taxonomy.tax",
            good_mags = "4500_assembly_analysis/good_mags.txt",
            proteom = "4500_assembly_analysis/proteomics/all_proteoms.faa",
            genome = "4500_assembly_analysis/genomics/all_genomes.fna",
            CDSs = "4500_assembly_analysis/proteomics/all_CDSs.fna",
#            gffs = "4500_assembly_analysis/proteomics/all_gffs.gff",
   run :
       from numpy import logical_and
       from os.path import join as pjoin

       # cutoff based on https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/
       print(len(input))

       tax_files = input
       comp = params.completness
       cont = params.contamination
       mag_stats = [f.replace("full_taxonomy.tax","magstats.csv") for f in tax_files]

       genome_fold = pjoin(os.getcwd(), os.path.dirname(output.magstats), "genomics")
       proteoms_fold = pjoin(os.getcwd(), os.path.dirname(output.magstats), "proteomics")

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
       os.makedirs(pjoin(proteoms_fold, "gffs"), exist_ok = True)
       os.makedirs(pjoin(proteoms_fold, "CDSs"), exist_ok = True)

       for f in tqdm(mag_stats):
           fold = pjoin(os.getcwd(), os.path.dirname(f))
           mag_list = os.listdir(pjoin(fold, "clean_bins"))
           for b in tqdm(mag_list):
                os.symlink(pjoin(fold, "clean_bins", b, b + ".fna"), pjoin(genome_fold, "genomes", b + ".fna"))
                os.symlink(pjoin(fold, "clean_bins", b, b + ".gff"), pjoin(proteoms_fold, "gffs", b + ".gff"))
                os.symlink(pjoin(fold, "clean_bins", b, b + ".ffn"), pjoin(proteoms_fold, "CDSs", b + ".ffn"))
                os.symlink(pjoin(fold, "clean_bins", b, b + ".faa"), pjoin(proteoms_fold, "proteoms", b + ".faa"))
       cat_line = "find {fold} -type l -exec cat {{}} \; > {file}"

       print("cat-ing proteom")
       call(cat_line.format(fold = pjoin(proteoms_fold, "proteoms"), file = output.proteom), shell= True)
       print("cat-ing genome")
       call(cat_line.format(fold = pjoin(genome_fold, "genomes"), file = output.genome), shell= True)
       print("cat-ing CDSs")
       call(cat_line.format(fold = pjoin(proteoms_fold, "CDSs"), file = output.CDSs), shell= True)

rule pre_cluster_proteom:
    params : id = config['analyses']['pre_cluster_proteom']['id'],
             length_cut = config['analyses']['pre_cluster_proteom']['length_cut'],
             temp_folder = config['general']['temp_dir']
    input : proteom = "{path}/proteomics/all_proteoms.faa"
    output : clusters = "{path}/proteomics/precluster/representative_seq.faa",
             groups = "{path}/proteomics/precluster/clusters.tsv"
    threads : 16
    run :
        groups = output.groups
        clusters = output.clusters
        temp_fasta = pjoin(params.temp_folder, "cdhit_temp.faa")
#        temp_fasta = pjoin(wildcards.path, "proteomics", "cdhit_temp.faa")
        print("running cd-hit")
        call("cd-hit -i {input} -o {output} -c {id} -M 0 -T {threads} -d 0 -s {cut}".format(input = input.proteom, output = temp_fasta, id = params.id, cut=params.length_cut, threads= threads), shell=True)
        clstrs = {}
        record = []
        print("parsing clusters")
        
        counter = 0
        with open(temp_fasta + ".clstr") as handle:
            for l in tqdm(handle):
                if l.startswith(">"):
                    if record != []:
                        clstrs[rep]  = (counter,record)
                        record = []
                    counter += 1
                else :
                    dat = l.split()[2]
                    if l.split()[-1] == "*" : 
                        rep = dat[1:][:-3]
                    record += [dat[1:][:-3]]
        clstrs[rep]  = (counter,record)

        print("outputing clusters")
        with open(groups, "w") as outp:
            outp.writelines(["Cluster_" + str(rec[1][0]) + "\t" + "\t".join(rec[1][1]) + "\n" for i,rec in enumerate(clstrs.items())])
        print("outputing fasta")
        
        with open(clusters, "w") as outp :
            counter = 0
            with open(temp_fasta) as handle:
                for l in tqdm(handle):
                    if l.startswith(">") :
                        if counter > 0 :
                            outp.writelines(record)
                        counter += 1
                        record = [ ">Cluster_" + str(clstrs[l.split()[0][1:]][0]) + "\n"]
                    else :
                        record += [l]


rule self_diamond :
    input : clusters = "{path}/proteomics/precluster/representative_seq.faa",
    output : hits =  "{path}/proteomics/precluster/selfhits.diamond"
    threads : 16
    shell :"""
diamond makedb --db {input.clusters} --in {input.clusters}
diamond blastp --more-sensitive  -e0.001  -p {threads} -f 6 -q {input.clusters} --db {input.clusters} -o {output.hits}
"""

rule silix :
    input : hits =  "{path}/proteomics/precluster/selfhits.diamond",
            clusters = "{path}/proteomics/precluster/representative_seq.faa",
            groups = "{path}/proteomics/precluster/clusters.tsv"
    output : raw_clusters = "{path}/proteomics/silix_clusters/clusters.tsv",
             clusters = "{path}/proteomics/cogs.tsv",
#             net =  "{path}/proteomics/precluster/selfhits.net",
    threads : 1
    run :
        from tqdm import tqdm
        preclust_file = input.groups
        clust_file = output.raw_clusters
        out_file = output.clusters

        os.system(" ".join(["silix ", input.clusters, input.hits, ">" , clust_file]))

        with open(preclust_file) as handle:
            preclust = {l.split()[0] : l.split()[1:] for l in tqdm(handle)}

        with open(clust_file) as handle:
            clust = {ll : l.split()[0] for l in tqdm(handle) for ll in preclust[l.split()[1]] }

        cog_ids = set(clust.values())
        cogs = {c : [] for c in cog_ids}
        for cds, cog in clust.items() :
            cogs[cog] += [cds]

        with open(out_file, "w") as handle:
            handle.writelines(["COG_" + k + "\t" + "\t".join(v) + "\n" for k, v in cogs.items()])

rule mags2cogs:
    input : clusters = "{path}/proteomics/cogs.tsv",
    output : mag2cogs = "{path}/mags/mag2cogs.tsv",
    threads : 1
    run :
        with open(input.clusters) as handle:
            cogs = {r.split()[0] : r.split()[1:] for r in handle.readlines()}
        cogs2mags = { k : {"_".join(vv.split("_")[:-1]) for vv in v} for k, v in cogs.items()}
        mags = set.union(*cogs2mags.values())
        mags2cogs = {m : [] for m in mags}
        for k, v in cogs2mags.items():
            for vv in v:
                mags2cogs[vv] += [k]
        mags2cogs = {k: set(v) for k, v in mags2cogs.items()}
        with open(output.mag2cogs, "w") as handle:
            handle.writelines([k + "\t" + "\t".join([vv for vv in v]) + "\n" for k,v in mags2cogs.items() ])

rule magNet:
    input : mags2cogs = "{path}/mags/mag2cogs.tsv",
    output : pairs = "{path}/mags/mag_pairs.csv",
    threads : 5
    run :
        with open(input.mags2cogs) as handle:
            mags2cogs = {r.split()[0] : set(r.split()[1:]) for r in handle.readlines()}
        handle = open(output.pairs, "w")
        o = []
        handle.writelines(["query\tsubject\tprop_of_q_in_s\n"])
        for k1, v1 in tqdm(mags2cogs.items()):
            l1 = len(v1)
            if "unbinned" not in k1 and l1 > 5 :
                for k2,v2 in mags2cogs.items():
                    if k2 == k1:
                        break
                    l2 = len(v2)
                    if "unbinned" not in k2 and l2 > 5  :
                        l2 = len(v2)
                        inter = len(v1.intersection(v2))
                        if inter/l1 > 0.3:
                            o += ["{p}\t{q}\t{d}\n".format(p = k1, q = k2, d = inter/l1)]
                        if inter/l2 > 0.3:
                            o += ["{p}\t{q}\t{d}\n".format(p = k2, q = k1, d = inter/l2)]
                handle.writelines(o)
                o = []
        handle.close()

rule ani_otus:
    params : simil = 95, 
    input : pairs = "{path}/mags/fastani_pairs.csv",
            magstats = "{path}/magstats.csv",
    output : otus = "{path}/mags/ani_otus.csv"
    threads : 5
    run :
        import igraph
        from numpy import  logical_and,logical_or
        import pandas

        pairs = input.pairs
        otus = output.otus
        simil = params.simil
        magstats = pandas.read_csv(input.magstats, index_col=0)
        decents = magstats.index[logical_and( magstats.completeness > 40, magstats.contamination <5)]
        subs = set(magstats.index ).difference(decents)
        with open(pairs) as handle:
            mag_sim = {tuple(l.split()[0:2]) : float(l[:-1].split()[2]) for l in tqdm(handle) if l[0] != "q"}

        def make_ani_otus(sim_dict, simil) : 
            valid_pair = lambda k : mag_sim.get((k[1],k[0])) and mag_sim[k] > simil and mag_sim[(k[1],k[0])] > simil
            good_pairs = [k for k in tqdm(sim_dict) if valid_pair(k)]
            species_graph = igraph.Graph()
            vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
            rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
            species_graph.add_vertices(len(vertexDeict))
            species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs if k[0] != k[1] ])
            genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]
            return genome_clusters

        genome_clusters = make_ani_otus(mag_sim, simil)

        left_pairs = {k : v for k, v in mag_sim.items() if k[0] not in decents and k[1] in decents and v > 95}
        subs = {l[0] : (None,0) for l in left_pairs}
        for p,ani in tqdm(left_pairs.items()):
            if subs[p[0]][1] < ani:
                subs[p[0]] = (p[1], ani)

        for k, v in subs.items():   
            for g in genome_clusters:
                if v[0] in g :
                    g += [k]
        
        with open(otus, "w") as handle:
            handle.writelines([ "aniOTU_{id}\t".format(id = i) + ";".join(gs) + "\n" for i, gs in enumerate(genome_clusters)])


rule ani_otus2:
    params : simil = 95, 
    input : pairs = "{path}/mags/fastani_pairs.csv",
            magstats = "{path}/magstats.csv",
    output : otus = "{path}/mags/ani_otus_{simil}_{comps}.csv"
    threads : 5
    run :
        import igraph
        from numpy import  logical_and,logical_or
        import pandas

        pairs = input.pairs
        otus = output.otus
        simil = float(wildcards.simil)
        comps = float(wildcards.comps)
        magstats = pandas.read_csv(input.magstats, index_col=0)
        decents = magstats.index[logical_and( magstats.completeness > comps, magstats.contamination <5)]

        with open(pairs) as handle:
            mag_sim = {tuple(l.split()[0:2]) : float(l[:-1].split()[2]) for l in tqdm(handle) if l[0] != "q"}

        def make_ani_otus(sim_dict, simil) : 
            valid_pair = lambda k : mag_sim.get((k[1],k[0])) and mag_sim[k] > simil and mag_sim[(k[1],k[0])] > simil and k[0] in decents and k[1] in decents
            good_pairs = [k for k in tqdm(sim_dict) if valid_pair(k)]
            species_graph = igraph.Graph()
            vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
            rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
            species_graph.add_vertices(len(vertexDeict))
            species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs if k[0] != k[1] ])
            genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]
            return genome_clusters

        genome_clusters = make_ani_otus(mag_sim, simil)        
        with open(otus, "w") as handle:
            handle.writelines([ "aniOTU_{id}\t".format(id = i) + ";".join(gs) + "\n" for i, gs in enumerate(genome_clusters)])



rule functional_otus:
    params : simil = 0.70, cog_count = 10
    input : pairs = "{path}/mags/mag_pairs.csv",
            mags2cogs = "{path}/mags/mag2cogs.tsv",
    output : otus = "{path}/mags/functional_otus.csv"
    threads : 5
    run :
        import igraph
        pairs = input.pairs
        otus = output.otus
        simil = params.simil
        cog_count = params.cog_count
        mags2cogs_file = input.mags2cogs

        with open(mags2cogs_file) as handle:
            mags2cogs = {r.split()[0] : set(r.split()[1:]) for r in handle.readlines()}

        with open(pairs) as handle:
            mag_sim = {tuple(l.split()[0:2]) : float(l[:-1].split()[-1]) for l in tqdm(handle) if l[0] != "q" if float(l[:-1].split()[-1]) > simil}

        def make_funct_otus(sim_dict, simil, cog_count ) : 
            valid_pair = lambda k : mag_sim.get((k[1],k[0])) and len(mags2cogs[k[0]]) > cog_count and len(mags2cogs[k[1]]) > cog_count 
            good_pairs = [k for k in tqdm(sim_dict) if valid_pair(k)]
            species_graph = igraph.Graph()
            vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
            rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
            species_graph.add_vertices(len(vertexDeict))
            species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs if k[0] != k[1] ])
            genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]
            return genome_clusters

        genome_clusters = make_funct_otus(mag_sim, simil, cog_count )
        mag2otu = {vv : i for i,v in enumerate(genome_clusters) for vv in v}
        
        genome_clusters = [['fMAG:' + gg for gg in g] for g in genome_clusters]
        left_pairs = {k : v for k, v in mag_sim.items() if k[0] not in mag2otu}

        subs = {l[0] : (None,0) for l in left_pairs}
        for p,ani in tqdm(left_pairs.items()):
            if subs[p[0]][1] < ani and p[1] in mag2otu:
                subs[p[0]] = (p[1], ani)

        for k, v in subs.items():   
            if v[1]:
                genome_clusters[mag2otu[v[0]]] += ['fSUB:' + k]


        
        with open(otus, "w") as handle:
            handle.writelines([ "cogOTU_{id}\t".format(id = i) + ";".join(gs) + "\n" for i, gs in enumerate(genome_clusters)])




#        mean = lambda l1,l2 : (l1+l2)/2
#        len_test = lambda l1, l2 : ((l1 < l2*1.1) and (l1> l2*0.9)) and ((l2 < l1*1.1) and (l2> l1*0.9))
#        for k1, v1 in tqdm(mags2cogs_sub.items()):
#            for k2,v2 in mags2cogs_sub.items():
#                if k2 == k1 :
#                    break
#
#                l1 = len(v1)
#                l2 = len(v2)

        #         if len_test(l1,l2) and ((k2, k1) not in pairs) :
        #             dist = len(v1.intersection(v2))/mean(l1,l2)
        #             if dist > 0.9:
        #                 pairs.add((k1,k2))

        # dists = {(k1,k2) : len(v1.intersection(v2))/len(v1) for k1, v1 in tqdm(mags2cogs_sub.items()) for k1, v2 in mags2cogs_sub.items() }

        # members = dict()
        # clusters = []
        # for p in pairs:
        #     if p[1] not in members and p[0] not in members:
        #         clusters.append([p[0], p[1]])
        #         cluster_id = len(clusters) -1
        #         members[p[0]] = cluster_id
        #         members[p[1]] = cluster_id
        #     elif p[1] not in members:
        #         members[p[1]] = members[p[0]]
        #         clusters[members[p[0]]] += [p[1]]
        #     else :
        #         members[p[0]] = members[p[1]]
        #         clusters[members[p[1]]] += [p[0]]
        # clusters = set([tuple(sorted([k] + v[1])) for k, v in dat.items() if len(v[1]) > 1])
        # with open("../full_taxonomy.tax") as handle:
        #     tax = {l.split(",")[0] : ";".join([ll for ll in l[:-1].split(",")[1:] if ll != ""]) for l in handle.readlines()}

        # taxed_mabs = {f : set([tax[vv] for vv in v if vv in tax] )for f,v  in mab_clusts.items()}

        # taxed_mabs = {k : v for k,v in taxed_mabs.items() if v }
        # mab_sizes = {k :  numpy.mean([len(mags2cogs[cc]) for cc in c]) for k,c in mab_clusts.items()}

        # def count_annot(bin, annot):
        #     with open(pjoin("proteoms/", bin + ".faa")) as handle :
        #         annots = [l for l in handle if l.startswith(">")]
        #     return len([l for l in annots if annot in l]) / len(annots)

        # tags = ['Photosys', 'NAD', 'hypothetical']
        # fast_class = lambda cluster : {t : numpy.mean([count_annot(c, t) for c in cluster]) for t in tags}
        # name_based_stats = { k : fast_class(v) for k, v in tqdm(list(mab_clusts.items())) if k not in taxed_mabs}
        # mab_gcs = {k : stats.GC[list(v)].mean() for k, v in tqdm(mab_clusts.items()) if k not in taxed_mabs}
        # stats = pandas.read_csv("../magstats.csv", index_col=0)
        # for k,v in name_based_stats.items():
        #     v['GC'] = mab_gcs[k]
        #     v['COG_count'] = mab_sizes[k]
        #     v['count'] = len(mab_clusts[k])

        # pandas.DataFrame.from_dict(name_based_stats, orient="index").to_csv("non_taxed_stats.csv")

rule fastANI_dists:
    params : block_size = 1070
    input : taxonomy = "{path}/full_taxonomy.tax",
    output : pairs = "{path}/mags/fastani_pairs.csv",
             paired = "{path}/mags/.paired"
    threads : 20
    run :
        import tempfile
        genome_fold = pjoin(os.path.dirname(input.taxonomy), "genomics", "genomes")
        with open(input.taxonomy) as handle :
            mags = [l.split(",")[0] for l in handle][1:]

        mag_blocks = [mags[i:(i+params.block_size)] for i in list(range(0,len(mags), params.block_size))]
        with open(output.pairs, "a") as handle:
            handle.writelines(["query\tsubject\tani\tsize_q\tsize_s\n"])

        for i,bloc1 in enumerate(mag_blocks):
            b1_tfile = tempfile.NamedTemporaryFile().name

            with open(b1_tfile, "w") as handle:
                handle.writelines([pjoin(genome_fold,l + ".fna") +"\n" for l in bloc1])
            for j,bloc2 in enumerate(mag_blocks):
                    print("doing bloc {i} and {j}".format(i = i, j=j))
                    b2_tfile = tempfile.NamedTemporaryFile().name
                    with open(b2_tfile, "w") as handle:
                        handle.writelines([pjoin(genome_fold,l + ".fna")  +"\n" for l in bloc2])

                    out_tfile = tempfile.NamedTemporaryFile().name
                    call("fastANI --ql {b1} --rl {b2} -o {out} -t {threads} #2> /dev/null".format(b1 = b1_tfile, b2 = b2_tfile, out = out_tfile, threads = threads), shell = True)
                    with open(out_tfile) as handle:
                        new_dat = ["\t".join([ll.split("/")[-1].replace(".fna","") for ll in l.split()]) +"\n" for l in handle.readlines()]
                    with open(output.pairs, "a") as handle:
                        handle.writelines(new_dat)

                    os.remove(out_tfile)
                    os.remove(b2_tfile)
            os.remove(b1_tfile)
            os.system("touch " + output.paired)

rule fastaANI_genome_clsts:
    params : ani_cutoff = 95
    input : taxonomy = "{path}/full_taxonomy.tax",
            pairs = "{path}/mags/fastani_pairs.csv",
            magstats = "{path}/magstats",
    output :
    threads : 20
    run :
        import igraph

        with open(input.pairs) as handle:
            pairs = { (l.split()[0],l.split()[1]) : float(l.split()[3]) for l in tqdm(handle.readlines()) if not l.startswith("query") and float(l.split()[3]) > params.ani_cutoff }
        species_graph = igraph.Graph()
        vertexDeict = { v : i for i,v in enumerate(set([x for k in pairs.keys() for x in k]))}
        rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
        species_graph.add_vertices(len(vertexDeict))
        species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k, v in hanis.items() if k[0] != k[1] ])
        genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.WEAK)]

rule make_cogs:
    input : cogs = "{path}/proteomics/cogs.tsv",
            proteom = "{path}/proteomics/all_proteoms.faa",
    output :  cogs = "{path}/cogs/cogs_metadata.tsv"
    threads : 1
    run :
            import pandas
            md = pandas.read_csv(input.magstats, index_col = 0)
            proteom_fold = pjoin(os.path.dirname(input.proteom), "proteoms")
            gffs_fold = pjoin(os.path.dirname(input.proteom), "gffs")
            cog_folder = pjoin(os.path.dirname(output.cogs), "cogs")
            os.makedirs(pjoin(cog_folder), exist_ok = True)

            with open(input.clusters) as handle:
                cogs = {r.split()[0] : r.split()[1:] for r in handle.readlines()}

            def find_in_fasta(id, file) :
                with open(file) as handle:
                    for ll in handle:
                        if ll.split()[0] == ">" + id:
                            break
                    header = ll[1:-1]
                    seq = ""
                    for ll in handle:
                        if ll.startswith(">"):
                            break
                        else :
                            seq += ll[:-1]
                return (header, seq)

            def find_in_gff(id, file) :
                with open(file) as handle:
                    for ll in handle:
                        vals = ll.split("\t")
                        if ll[0] != "#" and vals[8].startswith("ID=" + id):
                            break
                return vals

            cog_md = {}
            consensus_prod = lambda prods : [(p, prods.count(p)) for p in set(prods) if p != "hypothetical protein"]
            for name, cog in tqdm([c for c in cogs.items() if len(c[1]) >1 ]):
                faas = []
                gffs = []
                cog_md[name] = {}
                for g in cog:
                    bin = "_".join(g.split("_")[:-1])
                    faas += [find_in_fasta(g, pjoin(proteom_fold,bin +".faa"))]
                    gffs += [find_in_gff(g, pjoin(gffs_fold,bin +".gff"))]
                comment = [dict(tuple(ff.split("=")) for ff in f[8][:-1].split(";")) for f in gffs]
                prods = [c['product'] for c in comment]
                eCs = [c.get('eC_number', "hypothetical protein") for c in comment]
                symb = [c.get('Name', "hypothetical protein").split("_")[0] for c in comment]
                concencus = consensus_prod(prods)
                ec_concencus = consensus_prod(eCs)
                symb_concencus = consensus_prod(symb)
                if len(concencus) == 0 :
                    cog_md[name]['annot'] = "hypothetical protein"
                    cog_md[name]['hyp_prop'] = 1
                    cog_md[name]['annot_cert'] = 0
                else :
                    concencus = max(concencus, key=lambda x : x[1])
                    cog_md[name]['annot'] = concencus[0]
                    cog_md[name]['hyp_prop'] = (prods.count('hypothetical protein'))/len(cog)
                    cog_md[name]['annot_cert'] = concencus[1]/len(prods)

                if len(symb_concencus) == 0 :
                    cog_md[name]['symb'] = "NA"
                    cog_md[name]['no_symb'] = 1
                    cog_md[name]['symb_cert'] = 0
                else :
                    symb_concencus = max(symb_concencus, key=lambda x : x[1])
                    cog_md[name]['symb'] = symb_concencus[0]
                    cog_md[name]['no_symb'] = (symb.count('hypothetical protein'))/len(cog)
                    cog_md[name]['symb_cert'] = symb_concencus[1]/len(symb)

                if len(ec_concencus) == 0 :
                    cog_md[name]['eC_numb'] = "NA"
                    cog_md[name]['no_ec'] = 1
                    cog_md[name]['ec_cert'] = 0
                else :
                    ec_concencus = max(ec_concencus, key=lambda x : x[1])
                    cog_md[name]['eC_numb'] = ec_concencus[0]
                    cog_md[name]['no_ec'] = (eCs.count('hypothetical protein'))/len(cog)
                    cog_md[name]['ec_cert'] = ec_concencus[1]/len(eCs)

                cog_md[name]['seq_count'] = len(cog)
                with open(pjoin(cog_folder, name + ".faa"), "w") as handle:
                    handle.writelines(sum([[f[0] + " " + name + "\n", f[1] + "\n" ] for f in faas],[]))
