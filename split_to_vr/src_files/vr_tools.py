#!/usr/bin/env python3

from Bio import AlignIO
from Bio import SeqIO
import sys
from datetime import datetime
import os
import shutil


def write_vr(msa, genename, vr_dict, verbose=False):
    """
        input:
            1) a dictionary of the variable region start and stop coordinates
            2) the name of the gene
            3) a MSA to be diced up
        output:
            one file for each variable region, with gaps (the '-' character) removed
    """
    runtime=datetime.now().strftime('%m_%d_%H')
    
    temp_dir="../raw_data/{0}_temp_vr_{1}".format(genename, runtime)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
        if verbose==True:
            print("creating "+temp_dir)
        
    final_dir="../processed_data/{0}_variable_regions_{1}".format(genename, runtime)
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
        if verbose==True:
            print("creating "+final_dir)
            
    print("writing to folder {0}".format(final_dir))   
    for v,w in vr_dict.items():
        smaller=msa[:, w[0]:w[1]]

        with open("{2}/temp_{0}_{1}.fna".format(v, genename, temp_dir), "w") as f:
            AlignIO.write(smaller, f, "fasta")
                
        #ungap=SeqIO.parse("{2}/temp_{0}_{1}.fna".format(v, genename, temp_dir), "fasta")

        # I can't figure out how to do this through Biopython
        # so I'm going to do it the dumb way:
        # parse the fasta file to a dictionary while removing the gap characters
        # then write the dictionary to a new fasta file
        dumb=dict()
        for u in SeqIO.parse("{2}/temp_{0}_{1}.fna".format(v, genename, temp_dir), "fasta"):
            dumb[u.id]=u.seq.ungap('-')
       
        with open("{2}/{0}_{1}.fna".format(genename, v, final_dir), "w") as f:
            for k,v in dumb.items():
                blocks=[v[i:i+60] for i in range(0, len(v), 60)]
                f.write(">{0}\n".format(k))
                for b in blocks:
                    f.write(str(b)+"\n")
    
    # delete the temp files
    if os.path.isdir(temp_dir):
        #os.remove(temp_dir)
        shutil.rmtree(temp_dir)
        if verbose==True:
            print("removing {}".format(temp_dir))
    else: 
        print("  could not find temp file {} to delete".format(temp_dir))
                    
                    
def id_to_species_16s():
    """
        the MSA file '4_13_16s_ecoli.aln' had all the 
        species names stripped out, and the alignment 
        sequences are only identified by the cryptic 
        id numbers. This dict maps the id numbers back 
        to the genus-species names as listed in 
        the 'checked_16s_wEcoli.fasta' file. 
        
        Just don't throw that file out. Ever.
    """
    species_dict=dict()     
    for ecol in SeqIO.parse("../../../16s_tcoffee/checked_16s_wEcoli.fasta", "fasta"):
        species_dict[ecol.id]=ecol.description.split(';')[-1].replace(' ','_')
    return species_dict

     
if __name__ == "__main__":
    print("\n")
    print("+----------------------------------------------+")
    print("|                                              |")
    print("|                  vr_tools                    |")
    print("|                                              |")
    print("+----------------------------------------------+")    
    print("\nHello from vr_tools!\n\nThe following functions are available in this file\n")
    print("* write_vr")
    print(write_vr.__doc__)
    print("* id_to_species_16s")
    print(id_to_species_16s.__doc__)
    
    

