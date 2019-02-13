from subprocess import call
from os.path import join as pjoin
import os
import shutil

rule clean_metabat:
    input :
        clusters = "{path}/binning/metabat/clusters.txt",
        assembly = "{path}/assembly.fna",
    output :
        folder = "{path}/binning/metabat/bins"
    run :
        import os
        from Bio import SeqIO
        from os.path import join as pjoin
        from tqdm import tqdm

        with open(input.clusters) as handle :
            clsts = { l.split()[0] : l.split()[1] for l in handle}
        bad_ids = set(clsts.values())
        mapy = { c  : i  for i,c in tqdm(enumerate(bad_ids)) if c != 0 }
        mapy['0'] = "unbinned"
        clsts = {k : mapy[v] for k, v in clsts.items()}
        vvs = list(clsts.values())
        zero = len(str(len(mapy)))
        seqs = {c : [] for c in clsts.values() }
        seq_count = {c : 0 for c in clsts.values()}
        seq_tot = {c : vvs.count(c) for c in tqdm(set(vvs))}
        seq_zeros = {c : len(str(vvs.count(c))) for c in tqdm(set(vvs))}
#        print(list(clsts.keys())[0:10])
        for s in tqdm(SeqIO.parse(input.assembly, "fasta")):
            nam = s.id.split()[0]
            if nam in clsts:
                c = clsts[nam]
                b_id = str(c).zfill(zero)
                seq_count[c] += 1
                pos = str(seq_count[c]).zfill(seq_zeros[c]) + "/" + str(seq_tot[c])
                s.id = "bin-" + b_id + ":" + pos
                s.description = ""
                seqs[c] += [s]
        if not os.path.exists(output.folder):
            os.makedirs(output.folder)
        for k, v in seqs.items():
            SeqIO.write(v, pjoin(output.folder, "bin-" + str(k).zfill(zero) + ".fasta"), "fasta")




rule metabat :
    input : assembly = "{path}/filtered_assembly.fna",
            mapping = "{path}/filtered_assembly/mapping/map_table.tsv"
    output : file = "{path}/binning/metabat/clusters.txt"
    run :
        metabat_str = "metabat2 --maxP {maxP} --minS {minS} -m {min_len}  -s {min_bin_size} -i  {ass} -o {output} -a {mapping}  --saveCls  --unbinned -t {threads} --noBinOut"
        call(metabat_str.format(**config['binning']['metabat'], ass = input.assembly, output = output.file, mapping = input.mapping, threads = threads), shell = True)

rule concoct :
    params : min_len = 1500,
             min_bin_size = 10000,
             maxP = 93,
             minS = 50
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv"
    output : file = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}",
             concoct_abundances = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}.tsv"
    shell : """
    columns=`head -n1 {input.mapping} | sed 's/\t/\n/g' | grep -n bam | grep -v var | cut -f1 -d":" | sed 's/$/,/' | tr -d '\n'`
    cut -f1,${columns%%,} -d$'\t' {input.mapping} > {output.concoct_abundances}

    """

rule maxbin :
    params : max_exec = "/home/moritz/share/MaxBin-2.2.5/run_MaxBin.pl"
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv",
            assembly = "{path}/{name}/{assembler}/{name}.fna"
    output : file = "{path}/{name}/{assembler}/mapping/maxbin/maxbin_{name}",
             maxbin_abunds = "{path}/{name}/{assembler}/mapping/maxbin/maxbin_{name}.tsv"
    shell : """
    mkdir `dirname  {output.maxbin_abunds}`/abunds
    columns=`head -n1 {input.assembly} | tr "\t" "\n" | grep -n bam | grep -v var | cut -f1 -d":" | sed 's/$/,/' | tr -d '\n'`
    cut -f1,${{columns%%,}} -d$'\t' {input.mapping} > {output.maxbin_abunds}

    for f in `head -n1 {output.maxbin_abunds} | tr '\t' '\n'| grep -v contigName`
    do n=`head -n1 {output.maxbin_abunds} | tr "\t" "\n" | grep -n $f | cut -f1 -d":"`
    cut -f1,$n bla | grep -v contig > `basedir {output.maxbin_abunds}`/abunds/$f.tsv
    done

    ls `basedir {output.maxbin_abunds}`/abunds/$f.tsv > {output.maxbin_abunds}.lst


    {params.max_exec} -contig {input.assembly} -out {output.file} -abund_list {output.maxbin_abunds}.lst -thread {threads}
    """
