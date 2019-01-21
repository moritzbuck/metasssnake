from subprocess import call
from os.path import join as pjoin
import os
import shutil
import pandas
import hashlib
import numpy

shell.prefix("module load bioinfo-tools bbmap samtools perl_modules BioPerl prokka; ")

metasssnake_path = pjoin(os.environ['HOME'], "repos/moritz/metasssnake/")
config['general']['temp_dir'] = os.environ['SNIC_TMP']

if not os.path.exists(config['gtdb']['download']['local']):
    call("wget " + config['gtdb']['download']['remote'] + " -O " + config['gtdb']['download']['local'], shell=True)
if not os.path.exists(config['gtdb']['download']['refseq_local']):
    call("wget " + config['gtdb']['download']['refseq_remote'] + " -O " + config['gtdb']['download']['refseq_local'], shell=True)
if not os.path.exists(config['gtdb']['download']['genbank_local']):
    call("wget " + config['gtdb']['download']['genbank_remote'] + " -O " + config['gtdb']['download']['genbank_local'], shell=True)
if not os.path.exists(config['gtdb']['download']['uba_local']):
    shutil.makedirs(config['gtdb']['download']['uba_local'])
    call("wget " + config['gtdb']['download']['uba'] + " -O " + config['gtdb']['download']['uba_local'] + "../temp.tar.gz" , shell=True)
    call("tar xzvf " + config['gtdb']['download']['uba_local'] + "../temp.tar.gz" , shell=True)
    call("mv " + config['gtdb']['download']['uba_local'] + "../UBA*.fsa " + config['gtdb']['download']['uba_local'] , shell=True)
    os.remove(config['gtdb']['download']['uba_local'] + "../temp.tar.gz")

def get_taxa(wildcards):
    metadata = pandas.read_csv(config['gtdb']['download']['local'], sep = '\t', index_col = 0, low_memory=False)
    gtdb_ids = list(metadata.loc[[wildcards.taxon in c for c in  metadata.gtdb_taxonomy]].index)

    ncbi_fct = lambda g : "{first}/{second}/{third}/{gtdb_id}/pfams.json".format(first = g.split('_')[-1][0:3],second = g.split('_')[-1][3:6], third = g.split('_')[-1][6:9], gtdb_id = g)
    uba_fct = lambda g : "UBA/{gtdb_id}/pfams.json".format(gtdb_id = g)

    return [ncbi_fct(f) if not "UBA" in f else uba_fct(f) for f in gtdb_ids]

def all_taxa(wildcards):
    levels = {
        "phylum" : 1,
        "class" : 2,
        "order" : 3,
        "family" : 4,
        "genus" : 5,
        "species" : 6}
    metadata = pandas.read_csv(config['gtdb']['download']['local'], sep = '\t', index_col = 0, low_memory=False)
    all_tax = set([l.split(";")[levels[wildcards.level]] for l in metadata.gtdb_taxonomy])
    return ["{path}/{level}/{taxon}/ani_mat.csv".format(path = wildcards.path, level = wildcards.path, taxon =tax) for tax in all_tax]

rule download:
    output : gbk = "{path}/{gtdb_id}/genome.gbk",
             cdss = "{path}/{gtdb_id}/cdss.fna",
             genome = "{path}/{gtdb_id}/genome.fna",
             proteom = "{path}/{gtdb_id}/proteom.faa",
             metadata = "{path}/{gtdb_id}/metadata.json"
    params : temp_folder = pjoin(config['general']['temp_dir'], "{gtdb_id}")
    threads : 1
    run :
        def md5Checksum(filePath):
            with open(filePath, 'rb') as fh:
                m = hashlib.md5()
                for l in fh:
                    m.update(l)
                return m.hexdigest()
        metadata = pandas.read_csv(config['gtdb']['download']['local'], sep = '\t', index_col = 0, low_memory=False).loc[wildcards.gtdb_id].to_dict()


        if wildcards.gtdb_id.startswith("UBA"):
            dl_folder = params.temp_folder
            if not os.path.exists(dl_folder):
                os.makedirs(dl_folder)

            shutil.copy(pjoin( config['gtdb']['download']['uba_local'], wildcards.gtdb_id + ".fsa"), dl_folder )
            exe_str = "prokka --outdir {out_dir}  --force --prefix {prefix} --locustag {prefix} --cpus {threads} {infile}"
            call(exe_str.format(out_dir = dl_folder, prefix = wildcards.gtdb_id, threads = threads, infile = pjoin(dl_folder, wildcards.gtdb_id + ".fsa")), shell = True)
            genomics =  wildcards.gtdb_id + ".fna"
            proteomics =  wildcards.gtdb_id + ".faa"
            gbks =  wildcards.gtdb_id + ".gbf"
            cdss =  wildcards.gtdb_id + ".ffn"

            shutil.copy(pjoin(dl_folder, genomics), output.genome)
            shutil.copy(pjoin(dl_folder, proteomics), output.proteom)
            shutil.copy(pjoin(dl_folder, gbks), output.gbk)
            shutil.copy(pjoin(dl_folder, cdss), output.cdss)
        else :
            tries = 0
            checked = False
            while(tries < config['gtdb']['download']['retries'] and not checked):

                ncbi_id = metadata['ncbi_genbank_assembly_accession'].split(".")[0]
                refseq = pandas.read_table(config['gtdb']['download']['refseq_local'], skiprows=1, index_col='gbrs_paired_asm', low_memory=False)
                genbank = pandas.read_table(config['gtdb']['download']['genbank_local'], skiprows=1, index_col=0, low_memory=False)
                genbank.index = [t.split(".")[0] for t in genbank.index]
                refseq.index = [t.split(".")[0] for t in refseq.index]

                if ncbi_id in refseq.index:
                    ncbi_data = refseq.loc[ncbi_id].to_dict()
                elif ncbi_id in genbank.index :
                    ncbi_data = genbank.loc[ncbi_id].to_dict()
                else :
                    call("touch " + output.genome, shell=True)
                    call("touch " + output.proteom, shell=True)
                    call("touch " + output.gbk, shell=True)
                    call("touch " + output.cdss, shell=True)
                    sys.exit()

                metadata.update(ncbi_data)

                dl_folder = params.temp_folder 
                if not os.path.exists(dl_folder):
                    os.makedirs(dl_folder)
 
                call("wget -r -nd " + metadata['ftp_path'] + " -P " + dl_folder + " 2> /dev/null ", shell=True)

                all_files = os.listdir(dl_folder)

                md5_file = "md5checksums.txt"
                all_files = set(all_files)
                all_files.remove(md5_file)
                all_files.remove("README.txt")
                all_files.remove("assembly_status.txt")

                with open(pjoin(dl_folder, md5_file)) as handle :
                    md5s = {s.split()[1].split("/")[1] : s.split()[0] for s in handle}

                checked = all([md5Checksum(pjoin(dl_folder,f)) == md5s[f] for f in all_files if md5s.get(f)])
                tries += 1

            assert checked , "could not download after {tries}".format(tries = tries)

            genomics = [f for f in all_files if f.endswith("_genomic.fna.gz") and not "_from_genomic" in f]
            assert len(genomics) == 1, "More/less than one gneomics file?"
            genomics = genomics[0]

            proteomics = [f for f in all_files if f.endswith("_protein.faa.gz")]
            if len(proteomics) == 0 :
                call("unpigz " + pjoin(dl_folder, genomics), shell=True)

                exe_str = "prokka --outdir {out_dir}  --force --prefix {prefix} --locustag {prefix} --cpus {threads} {infile}"
                call(exe_str.format(out_dir = dl_folder, prefix = wildcards.gtdb_id, threads = threads, infile = pjoin(dl_folder, genomics[:-3])), shell = True)
                genomics =  wildcards.gtdb_id + ".fna"
                proteomics =  wildcards.gtdb_id + ".faa"
                gbks =  wildcards.gtdb_id + ".gbf"
                cdss =  wildcards.gtdb_id + ".ffn"

                shutil.copy(pjoin(dl_folder, genomics), output.genome)
                shutil.copy(pjoin(dl_folder, proteomics), output.proteom)
                shutil.copy(pjoin(dl_folder, gbks), output.gbk)
                shutil.copy(pjoin(dl_folder, cdss), output.cdss)
            else :
                assert len(proteomics) == 1, "More/less than one proteins seq file?"
                proteomics = proteomics[0]

                gbks = [f for f in all_files if f.endswith("_genomic.gbff.gz")]
                assert len(gbks) == 1, "More/less than one gneomics file?"
                gbks = gbks[0]

                cdss = [f for f in all_files if f.endswith("_cds_from_genomic.fna.gz")]
                assert len(cdss) == 1, "More/less than one gneomics file?"
                cdss = cdss[0]

                umzip = "unpigz -c {file} > {outfile}"
                call( umzip.format(file = pjoin(dl_folder, genomics), outfile = output.genome), shell=True)
                call( umzip.format(file = pjoin(dl_folder, proteomics), outfile = output.proteom), shell=True)
                call( umzip.format(file = pjoin(dl_folder, gbks), outfile = output.gbk), shell=True)
                call( umzip.format(file = pjoin(dl_folder, cdss), outfile = output.cdss), shell=True)

        metadata['genome'] = output.genome
        metadata['proteome'] = output.proteom
        metadata['gbk'] = output.gbk
        metadata['cdss'] = output.cdss

        metadata = { k : int(v) if type(v) == numpy.int64 else v for k,v in metadata.items()}

        with open(output.metadata, "w") as handle:
            json.dump(metadata, handle)

        shutil.rmtree(dl_folder)

rule hammer :
    input : proteom = "{path}/{gtdb_id}/proteom.faa"
    output : pfams = "{path}/{gtdb_id}/pfams.json", raw_out = "{path}/{gtdb_id}/raw_pfams.tblout"
    run :
        import pandas

        call("hmmsearch --noali --tblout {raw_out}  {db} {input}  > /dev/null".format(raw_out = output.raw_out, db = config['gtdb']['hmmer']['pfams_db'], input = input.proteom), shell = True)
        domtblout_head = ["target_name" , "target_accession" , "query_name" , "query_accession" , "E-value","score-sequence" , "bias-sequence" , "bdE-value","score-best-domain" , "bias--best-domain" , "exp" , "reg" , "clu" , "ov" , "env" , "dom" , "rep" , "inc" , "description_of_target"]
        data = pandas.read_csv(output.raw_out, delim_whitespace=True, comment="#", names=domtblout_head[:-1], usecols=range(0,18))
        data = data.loc[data['bdE-value'] < 10e-6]
        pfams_tuples = { ( b['query_accession'], b['target_name'] ) for a,b in data.iterrows()}
        pfam_dict = {pfam : [b for a,b in pfams_tuples if a == pfam] for pfam in {a for a,b in pfams_tuples}}

        with open(output.pfams, "w") as handle:
            json.dump(pfam_dict, handle)

rule process_taxon:
    input : get_taxa
    output : pfam_matrix = "{path}/{taxon}/pfam_mat.csv",
             ani_matrix = "{path}/{taxon}/ani_mat.csv"
    run :
        call("touch {file}".format(file = output.pfam_matrix), shell = True)
        call("touch {file}".format(file = output.ani_matrix), shell = True)

rule all_taxa:
    input : all_taxa
    output : touch("{path}/{level}/all_computed.flag")
