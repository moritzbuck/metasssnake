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



rule clean_metabat:
    input : "{path}/{name}/{assembler}/mapping/metabat/metabat_{name}"
    output : "{path}/{name}/{assembler}/MAGs/metabat_{name}-unbinned.fa"
    run :
        import os
        from Bio import SeqIO
        from os.path import join as pjoin
        from tqdm import tqdm

        ipath = pjoin(os.path.dirname(input[0]))
        opath = output[0]
        os.makedirs(opath)
        for f in tqdm(os.listdir(ipath)):
            if f[-3:] == ".fa":
                with open(pjoin(ipath, f)) as handle:
                    seqs = [s for s in SeqIO.parse(handle, "fasta")]
                zeros = len(str(len(seqs)))
                bin_name = f[:-3]
                for i,s in enumerate(seqs):
                    s.id = bin_name.replace(".","-") + "-" + str(i+1).zfill(zeros)
                    s.description = ""
                with open(pjoin(opath, f[:-3].replace(".","-")+".fa"), "w") as handle:
                    SeqIO.write(seqs, handle, "fasta")


rule metabat :
    params : min_len = 1500,
             min_bin_size = 10000,
             maxP = 93,
             minS = 50
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv"
    output : file = "{path}/{name}/{assembler}/mapping/metabat/metabat_{name}"
    threads : THREADS
    shell : """
    metabat2 --maxP {params.maxP} --minS {params.minS} -m {params.min_len}  -s {params.min_bin_size} -i  {wildcards.path}/{wildcards.name}/{wildcards.assembler}/{wildcards.name}.fna -o {output.file} -a {input.mapping}  --saveCls  --unbinned -t {threads}
    """

rule concoct :
    params : min_len = 1500,
             min_bin_size = 10000,
             maxP = 93,
             minS = 50
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv"
    output : file = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}",
             concoct_abundances = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}.tsv"
    threads : THREADS
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
    threads : THREADS
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
