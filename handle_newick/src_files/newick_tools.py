#!/usr/bin/env python3

from Bio import Phylo
from Bio import SeqIO
import math
import statistics
import sys
from datetime import datetime
import os
import ast
import collections

def update_progress(progress, what_doing):
    """
        update_progress() : Displays or updates a console progress bar
        Accepts a float between 0 and 1. Any int will be converted to a float.
        A value under 0 represents a 'halt'.
        A value at 1 or bigger represents 100%

        Shamelessly copied and pasted from stackoverflow
        because it's very cool. Written by Brian Khuu
        https://stackoverflow.com/users/2254146/brian-khuu
    """
    barLength = 30 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\r{3}: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100), status, what_doing)
    sys.stdout.write(text)
    sys.stdout.flush()

def id_to_species_16sOLD():
    """
        replaced with inhouse_to_species()
        
        the MSA file 'thesis/16s_tcoffee/4_13_16s_ecoli.aln' had all the 
        species names stripped out, and the alignment 
        sequences are only identified by the cryptic 
        id numbers. This dict maps the id numbers back 
        to the genus-species names as listed in 
        the reference file, which can either be bladder or vaginal. 
        
        made a copy of checked_16s_wEcoli.fasta 
        and put it into the thesis_git/resources/ directory
    """
    # bladder species
    path_to_file="../../resources/bladder_16s_silva.fasta"
    
    # vaginal species
    #path_to_file="../../resources/checked_16s_wEcoli.fasta" 
    
    species_dict=dict()     
    for ecol in SeqIO.parse(path_to_file, "fasta"):
        #species_dict[ecol.id]=ecol.description.split(';')[-1].replace(' ','_')
        species_dict[ecol.id]=ecol.description.split(';')[-1]
    return species_dict
    
def inhouse_to_species(path_to_file):
    """
        The record description of a fasta file is lost along the pipeline. 
        inhouse_to_species() returns a dictionary that maps accession numbers 
        to the taxonomy that was originally with the record. 

        Some map files used so far are listed below.

            TKW dataset with Wellcome-Trust IDs     /resources/wellcome_to_species.csv
            TKW dataset of downloaded species       /resources/all_16s_bladder_5_8.csv
            
            older files not used anymore:
                                                    /resources/master_strain_list_05_07.csv                                                   

        You can pass a file of your own. The csv must be formatted like this:

            16933_8_1,Corynebacterium coyleae
            16933_8_2,Gardnerella vaginalis
            ...
    """
    
    species_dict=dict()   
    with open(path_to_file, 'r') as f:
        for line in f:
            (k,v)=line.split(',')
            species_dict[k]=v.strip()
    return species_dict
    
def vr_per_gene(route):
    """
        need a way to count how many variable regions 
        are listed for each gene. Using a copy of the demarcated 
        variable regions file from 
        sliding_window/processed_data/04_12_18_29_vr_width_output/04_12_automap_vr_cut463.csv
        
        returns dictionary of {gene name : number of variable regions}
    """
    
    with open(route, 'r') as s:
        # read from a generator into dictionary
        # using 'ast' module to convert a string written like a list to an actual list
        # thank you, fraxel from stackoverflow           
        get_coords={l.split('\t')[0].strip():len(ast.literal_eval(l.split('\t')[1].strip())) for l in s}
        #get_coords=[l.split('\t')[0].strip() for l in s]
        
    return get_coords
    
def vr_leaf_dict(vr_values, path_to_folder, prefix_of_file):
    """
        returns dictionary of {vr:Phylo object of newick guide tree}
        
        hard coded file names! change that!
    """
    tree_dict=dict()
    for x in range(1, vr_values+1, 1):

        #tree_dict["v"+str(x)]=Phylo.read("{0}/{1}_{2}.dnd".format(path_to_folder, prefix_of_file, "v"+str(x)), "newick")
        tree_dict["v"+str(x)]=Phylo.read("{0}/{1}_{2}_nostrains.dnd".format(path_to_folder, prefix_of_file, "v"+str(x)), "newick")

    return tree_dict
    
def tree_success(vrs_in_gene, cutoff, new_dir, mg_or_16s, routing, file_prefix):
    """
        inputs
            vrs_in_gene    a dictionary of vr:newick guide tree
            cutoff         cutoff value of branch distance
            new_dir        path to a directory to hold the file
            mg_or_16s      running functional or ribosomal gene?
            routing        path to vr .dnd files
            file_prefix    gene or ribo file prefix?
            
        output
            writes the percent success of each vr
            of a marker gene to correctly identify 
            the species greater than a cutoff to compared
            to the full guide tree
    """
    #translation=id_to_species_16s()
    #mg_translation=inhouse_to_species() 
    
    ribo_translation=inhouse_to_species("../../resources/all_16s_bladder_5_8.csv")
    mg_translation=inhouse_to_species("../../resources/wellcome_to_species.csv")
    silva_to_tkw=inhouse_to_species("../../resources/silva_tkw_species_translation.csv")
    
    if mg_or_16s:
        pf="functional"
    else:
        pf="ribo"
        
    # write header: gene  vr  #identified  #missed  #strains  #total  %success
    with open("{0}/{1}_success_results.txt".format(new_dir,pf), "a") as q:
        q.write("gene\tvr\tidentified\tmissed\ttotal\tpercent_success\n")
        
    # iterate through all the variable region guide trees
    for k,v in vrs_in_gene.items():

        if mg_or_16s==True:
            # use the gene names as file prefixes
            tree_dict=vr_leaf_dict(v, routing, k)
        else:
            # use the 16s bladder names as prefix
            tree_dict=vr_leaf_dict(v, routing, file_prefix)   

        for m,n in tree_dict.items():
         
            twigs=[i for i in n.find_clades(terminal=True) if (i.branch_length>=cutoff)]
            #strains=[i for i in n.find_clades(terminal=True) if (i.branch_length<.001)]
            check=[i for i in n.find_clades(terminal=True)]
            
            percent_accurate=len(twigs)/len(check)
            missed=len(check)-len(twigs)

            if mg_or_16s:
                twig_names=[mg_translation[tn.name[:-6]] for tn in twigs] 
                #strain_names=set([mg_translation[sn.name[:-6]] for sn in strains])
                check_names=set([mg_translation[cn.name[:-6]] for cn in check])  
            else:
                twig_names=[ribo_translation[tn.name] for tn in twigs]
                #strain_names=set([ribo_translation[sn.name] for sn in strains])
                check_names=set([ribo_translation[cn.name] for cn in check])  
             
            #print("\n*** building {}".format(m))
            with open("{0}/{1}_{2}_matrix.csv".format(new_dir,pf,m), "a") as v: 
                v.write("species,{}\n".format(m))
                for c in sorted(check_names):
                    if c in twig_names:
                        v.write("{0},1\n".format(silva_to_tkw[c]))
                        
                    else:
                        v.write("{0},0\n".format(silva_to_tkw[c]))

            
            #with open("{0}/{1}_success_results.txt".format(new_dir,pf), "a") as q: 
                # write to file: gene  vr  #identified  #missed  #strains  #total  %success
            #    q.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k, m, len(twigs), missed, len(check), percent_accurate))

def remove_strainsOLD(path_to_dir, filename, output_dir, mg_16s, verbose=False):
    """        
        the old version of this function.
        Generates a strainless file by building a new dictionary every time the function is called. 
        I thought this was a bad idea because python dictionaries retrieve key:value pairs in different order than when entered
    """

    mg_translation=inhouse_to_species("../../resources/master_strain_list_05_07.csv") 
    translation=id_to_species_16s()
    gene=filename[:-4]

    # read in the tree
    #temp_tree=Phylo.read("../raw_data/vr_w_aln4_26/dnd/COG0090_v3.dnd", "newick")
    temp_tree=Phylo.read("{0}/{1}".format(path_to_dir, filename), "newick")

    if mg_16s:
        # need to make another dictionary that maps the in-house isolate name 
        # with the wellcome trust addition of 6 characters. What a pain in the ass.
        all_isolates=[i.name for i in temp_tree.find_clades(terminal=True)]
        wellcome_name={s[:-6]:s for s in all_isolates}

        invert=[(v,k) for k,v in mg_translation.items()]
    else:
        invert=[(v,k) for k,v in translation.items()]

    # from the defaultdict docs:
    #    When each key is encountered for the first time, 
    #    it is not already in the mapping; so an entry is 
    #    automatically created using the default_factory 
    #    function which returns an empty list. The list.append() 
    #    operation then attaches the value to the new list. When 
    #    keys are encountered again, the look-up proceeds normally 
    #    (returning the list for that key) and the list.append() 
    #    operation adds another value to the list. 
    #
    # so that's nice
    strain_catalogue=collections.defaultdict(list)
    for k,v in invert:
        strain_catalogue[k].append(v)

    # this iterates through the entire mg_translation dictionary   
    for k,v in strain_catalogue.items():
        if mg_16s:
            batch=["{}_{}".format(mg_translation[b],b) for b in v]
        else:
            batch=["{}_{}".format(translation[b],b) for b in v]
        if verbose:
            print(batch)
        with open("../processed_files/removed_strains_{0}_output.txt".format(gene,output_dir), "a") as f:
            f.write("\nkeeping {}\n".format(batch[0]))
            # only remove those isolates that were identified as strains
            if len(v)>1:       
                for x in v[1:]:     
                    if mg_16s:
                        f.write("\tremoving {}_{}\n".format(mg_translation[x], x))
                        temp_tree.prune(target=wellcome_name[x])
                    else:
                        f.write("\tremoving {}_{}\n".format(translation[x], x))
                        temp_tree.prune(target=x)

    Phylo.write(temp_tree, "{1}/{0}_nostrains.dnd".format(gene, output_dir), "newick")
    
def remove_strains(path_to_dir, filename, output_dir, mg_16s, verbose=False):
    """        
        Removes the strains from a newick tree, based on reference files.
        
        input
            path_to_dir  path to directory with newick files
            filename     the filename inside that directory
            output_dir   path to the folder to write the results
            mg_16s       flag that indicates functional (True) or ribosomal sequences (False)
            
        output
            a new newick file with the strains removed, written to /processed_data
            a log of the strains removed, written to /processed_files
    """

    mg_translation=inhouse_to_species("../../resources/wellcome_to_species.csv") 
    ribo_translation=inhouse_to_species("../../resources/all_16s_bladder_5_8.csv")
    
    if mg_16s:
        strain_translation=inhouse_to_species("../../resources/master_strain_list_05_07.csv") 
    else:
        strain_translation=inhouse_to_species("../../resources/master_16s_bladder_strains_05_08.csv")

    gene=filename[:-4]
    
    strain_list=list(strain_translation.keys())
    runtime=datetime.now().strftime('%m_%d_%H_%M')

    # read in the tree
    temp_tree=Phylo.read("{0}/{1}".format(path_to_dir, filename), "newick")
    tree_ids=[i.name for i in temp_tree.find_clades(terminal=True)]

    with open("../processed_files/removed_strains_{0}_output{1}.txt".format(gene,runtime), "a") as f:   
        for x in tree_ids:
            # this is a sloppy hack, fix it. Stop farting around with the different ID names
            if mg_16s:
                y=x[:-6]
            else:
                y=x
                
            if y in strain_list:
                # only remove those isolates that were identified as strains
                if mg_16s:
                    f.write("\tremoving {} {}\n".format(mg_translation[x[:-6]], x[:-6]))
                    temp_tree.prune(target=x)
                else:
                    f.write("\tremoving {} {}\n".format(ribo_translation[x], x))
                    temp_tree.prune(target=x)
            else:
                if mg_16s:
                    f.write("keeping {} {}\n".format(mg_translation[x[:-6]], x[:-6]))
                else:
                    f.write("keeping {} {}\n".format(ribo_translation[x], x))

    Phylo.write(temp_tree, "{1}/{0}_nostrains.dnd".format(gene, output_dir), "newick")
    
def chop_branches(coord_path, dendro_path, verbose=False):
    """
        input
            the path to the file that lists how the dendrogram 
                is to be divided into subgroups
            the path to the full dendrogram generated by tcoffee
            
        output
            writes a new file with named internal nodes where 
                the new subgroups are descended
            returns dictionary of leaves in each subgroup
    """
    
    runtime=datetime.now().strftime('%m_%d_%H_%M')
    
    with open(coord_path, 'r') as s:
        # read from a generator into dictionary
        # using 'ast' module to convert a string written like a list to an actual list
        # thank you, fraxel from stackoverflow        
        get_roots={l.split(' ')[0].strip():ast.literal_eval(l.split(' ')[1].strip()) for l in s} 
            
    full_dendrogram=Phylo.read(dendro_path, "newick")

    # locates the common ancestor node in the full tree
    #alpha=full_dendrogram.common_ancestor(['17957_1_66_01215','16933_8_10_00286'])
    #beta=full_dendrogram.common_ancestor(['16933_8_11_01078','21837_8_92_00097'])
    #gamma=full_dendrogram.common_ancestor(['16933_8_6_01788','17957_1_46_04458'])

    common_node=dict()
    for k,v in get_roots.items():
        if verbose:
            print("common_node[{0}]=full_dendrogram.common_ancestor({1})".format(k,v))
        common_node[k]=full_dendrogram.common_ancestor(v)
        

    subdivided_tree=dict()
    for k,v in common_node.items():
        # names those nodes in case it's needed later
        v.name=k
        # returns all the leaves in the named clade 'gamma'
        all_leaves=[i.name for i in full_dendrogram.find_elements(terminal=True) if v.is_parent_of(i)]
        subdivided_tree[k]=all_leaves
        
    # this seems to write just fine, but UGENE throws an error about
    # Subtask {Load document: 'cladenames_201.nwk'} is failed: Error parsing weight: alpha
    # looked at the file and I can't see anything weird
    # ...   :0.10178):0.00699)alpha:0.00017,(((((((((((((16933_8_11_01078:0.0 ... <- written file
    # ...   :0.10178):0.00699):0.00017,(((((((((((((16933_8_11_010778:0.0 ... <- original file
    # not that it matters. I'm not using UGENE to manipulate the trees

    # writes a new tree that includes the names
    Phylo.write(full_dendrogram, "../processed_files/cog201_named_clades{0}.dnd".format(runtime), "newick")
    
    return(subdivided_tree)
    
if __name__== "__main__":

    print("\n")
    print("+----------------------------------------------+")
    print("|                                              |")
    print("|               newick_tools                   |")
    print("|                                              |")
    print("+----------------------------------------------+")    
    print("\nHello from newick_tools!\n\nThe following functions are available in this file\n")
    print("* update_progress(progress, what_doing)")
    print(update_progress.__doc__)
    print("* id_to_species_16s()")
    print(id_to_species_16s.__doc__)
    print("* inhouse_to_species(path_to_file='../../resources/wellcome_to_species.csv')")
    print(inhouse_to_species.__doc__)
    print("* vr_per_gene(route)")
    print(vr_per_gene.__doc__)
    print("* vr_leaf_dict(vr_values, path_to_folder, prefix_of_file)")
    print(vr_leaf_dict.__doc__)
    print("* tree_success(vrs_in_gene, cutoff, new_dir, mg_or_16s, routing, file_prefix)")
    print(tree_success.__doc__)
    print("* remove_strains(path_to_dir, filename, output_dir, mg_16s, verbose=False)")
    print(remove_strains.__doc__)
    print("* chop_branches(coord_path, dendro_path, verbose=False)")
    print(chop_branches.__doc__)


