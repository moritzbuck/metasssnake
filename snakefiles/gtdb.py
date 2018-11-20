from subprocess import call
from os.path import join as pjoin
import os
import shutil

metasssnake_path = "/home/moritzbuck/repos/moritz/metasssnake/"
configfile : pjoin(metasssnake_path, "snakefiles" , "params.json")


rule download:
    output : gbk = "{path}/{gtdb_id}/genome.gbk", cdss = "{path}/{gtdb_id}/cdss.fna", genome = "{path}/{gtdb_id}/genome.fna", proteom = "{path}/{gtdb_id}/proteom.faa", metadata = "{path}/{gtdb_id}/metadata.json"
    params : temp_folder = pjoin(config['general']['temp_dir'], "{gtdb_id}")
    threads : 1
    run :
        import pandas
        import hashlib
        import numpy

        def md5Checksum(filePath):
            with open(filePath, 'rb') as fh:
                m = hashlib.md5()
                for l in fh:
                    m.update(l)
                return m.hexdigest()

        tries = 0
        checked = False
        while(tries < config['gtdb']['download']['retries'] and not checked):
            if not os.path.exists(config['gtdb']['download']['local']):
                call("wget " + config['gtdb']['download']['remote'] + " -O " + config['gtdb']['download']['local'], shell=True)
            if not os.path.exists(config['gtdb']['download']['refseq_local']):
                call("wget " + config['gtdb']['download']['refseq_remote'] + " -O " + config['gtdb']['download']['local'], shell=True)
            if not os.path.exists(config['gtdb']['download']['genbank_local']):
                call("wget " + config['gtdb']['download']['genbank_remote'] + " -O " + config['gtdb']['download']['local'], shell=True)

            metadata = pandas.read_csv(config['gtdb']['download']['local'], sep = '\t', index_col = 0, low_memory=False).loc[wildcards.gtdb_id].to_dict()
            ncbi_id = metadata['ncbi_genbank_assembly_accession']
            refseq = pandas.read_table(config['gtdb']['download']['refseq_local'], skiprows=1, index_col='gbrs_paired_asm', low_memory=False)
            genbank = pandas.read_table(config['gtdb']['download']['genbank_local'], skiprows=1, index_col=0, low_memory=False)
            if ncbi_id in refseq.index:
                ncbi_data = refseq.loc[ncbi_id].to_dict()
            else :
                ncbi_data = genbank.loc[ncbi_id].to_dict()
            metadata.update(ncbi_data)

            dl_folder = pjoin(config['general']['temp_dir'], wildcards.gtdb_id)
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
