import os
import shutil
import tempfile

from get_seqs_functions import get_ref_by_gene

def get_ref_by_gene_resolve(fout = None, no_bed = False, directory = None, **kwargs):
    
    with tempfile.TemporaryDirectory() as tmpdir:
        if not directory:
            directory = os.getcwd()
        get_ref_by_gene(**kwargs, directory = directory)
        
        ## merge fasta files if required
        if fout is not None:
            seqs_out = {}
            from fasta_manip import fasta_to_dict, dict_to_fasta
            for fname in os.listdir(tmpdir):
                if fname[-6:] == ".fasta":
                    seqs_out = {**seqs_out, **fasta_to_dict(os.path.join(tmpdir, fname))}
            dict_to_fasta(seqs_out, fout)
        
        ## mv bed directory to output directory if no_bed not raised
        if not no_bed:
            out_dir = os.path.join(directory, "bed")
            os.makedirs(out_dir, exist_ok = True)
            for fname in os.listdir(os.path.join(tmpdir, "bed")):
                shutil.move(os.path.join(tmpdir, bed, fname), os.path.join(out_dir, fname))
    return
