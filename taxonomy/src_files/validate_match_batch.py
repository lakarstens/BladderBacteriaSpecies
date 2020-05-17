#!/usr/bin/env python3

import sys
# since I want to import a file from another folder, 
# I'm going to use sys.path.insert to tell 
# python where to look at runtime. As seen on stackoverflow.
sys.path.insert(0, "../../handle_newick/src_files")
import newick_tools
import argparse
import os
from datetime import datetime
import collections

#--------------------------
#       argparse!
#--------------------------

parsing=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""

    *** updating this script ***
    
    This script compares the results of the taxonomy assigned by BLCA to the known taxonomy held in a map file.
    It will also check for known synonyms of the result taxonomy, but that list is incomplete 
    until I find a better way to keep it updated. Misspellings that don't map to known synonyms are assigned 
    as a mismatch. I'll fix that later, but know that unless you're typing your own taxonomy into fasta records, 
    the chance of misspellings are low.
    
    Any genus or species taxon labled as "unidentified" or "uncultured" is designated as "unavailable". 
    "Unclassified" is reserved for when the classification algo thinks there's no match for the record. 
    "Unavailable" just means there's no information for that record, and I'm choosing to put it in 
    the True Non-Match box. In a test case environment, these records don't provide any way of determining 
    if the records match or not, so I think I have to throw them out. Databases like greengenes and silva 
    will suffer, because there's frequently no species taxon provided.
    
    Synonyms
    ------------------------------
    Prokaryotic synonyms were downloaded from the Prokaryotic Nomenclature Up-to-Date website run by the 
    DSMZ-German Collection of Microorganisms and Cell Cultures, GmbH. If the synonym names are already the same, 
    it's because the DSMZ included information down to the strain level. Strains are out of scope, so entries 
    like 'Enterobacter cloacae cloacae' and 'Enterobacter cloacae dissolvens' are considered synonyms of 'Enterobacter cloacae'
    
    Use
    ------------------------------
    Input is the formatted csv file as written by the python script 'blcaout_to_csv.py'
    
        id,taxon,rank,confidence
        CP011325.22226.23780,Bacteria,domain,100.0
    
    the --map_file option takes a path to a file that can be 
    used to map accession numbers to taxonomy by newick_tools.inhouse_to_species(). 
    See below for the list of some existing map files. 
    
        KTW sequencing data                     master_wellcome_list_2019-05-07
        Type strains of bladder KTW dataset     type_strains_map_2019-11-26
    
    run the script like this:
        
        >python validate_match.py -i path/to/formatted/folder -m path/to/map/file
        
    *** disabled for now ***
    and if you want to rename the file, write out the full or relative path but leave off the file type. 
    "_validated_<date>.csv" will be added to the end of the pathname.
    
        >python validate_match.py -i path/to/formatted/folder -m path/to/map/file -o ../relative/path/and/file_prefix
        
""", epilog="This script depends on the following *newick_tools* functions:\n  inhouse_to_species()\n{}".format(newick_tools.inhouse_to_species.__doc__))

parsing.add_argument("--infolder", "-i", required=True, help="The full or relative path to the folder.")
parsing.add_argument("--map_file", "-m", required=True, help="what accession:species map file will be used?")
parsing.add_argument("--outfolder", "-o", action="store_true", help="the infile name is probably really awkward, so you might want to rename it. But you must provide the FULL or RELATIVE PATH")
parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
parsing.add_argument("--debug", action="store_true", help="excruciating detail, sometimes written to a file")
  
args=parsing.parse_args()

#--------------------------
#   get to work
#--------------------------

# build the synonym dictionary
synonym_file="../../resources/current_synonyms_2020-01-05_30.csv"
synonym_dict=dict()
with open(synonym_file, 'r') as s:
    getlines=s.readlines()
    for lines in getlines:
        lineitems=lines.strip().split(',')
        synonym_dict[lineitems[0]]=lineitems[1]
    
runtime=datetime.now().strftime('%Y-%m-%d_%M')

get_files=os.listdir(args.infolder)

if args.debug:
    print(get_files)

if args.outfolder:
    outfile="{}_validated_{}.csv".format(args.outfile, runtime)
else:
    outfolder="../processed_files/validated_outfiles_{}".format(runtime)
    
#if args.outfile:
#    outfile="{}_validated_{}.csv".format(args.outfile, runtime)
#else:
#    outfile="../processed_files/think_of_a_better_name_validated_{}_{}.csv".format(os.path.basename(infolder).split('.')[0], runtime)


if not os.path.exists(outfolder):
    os.makedirs(outfolder)
    


# a dictionary of the dataset accession number and known taxonomy 
translate_dataset=newick_tools.inhouse_to_species(args.map_file)
if args.debug:
    print(translate_dataset)

# these are the problem words in classification
#bad_words=["unavailable", "uncultured", "metagenome", "unidentified"]
# on second thought, these are the bad words
bad_words=["uncultured", "unidentified", "unavailable"]


if args.verbose:
    print("reading in {}".format(args.infolder))
    
for gf in get_files:
    outfile="{2}/validated_{1}_{0}.csv".format(runtime, gf.split("/")[-1].split(".")[0].split("-")[0], outfolder)
    
    # make dictionary of blca results
    results_dict=dict()
    with open("{}/{}".format(args.infolder, gf), 'r') as f:
        filelines=f.readlines()
        # skip the first line, it's the colnames
        for lines in filelines[1:]:
            lineitems=lines.strip().split(',')
            if lineitems[2]=='genus':
                # gg sometimes has blank spaces in this taxon
                # empty strings evaluate to 'false'
                if not lineitems[1]:
                    genus="unavailable"
                # silva includes a number sometimes -> 'Corynebacterium 1'
                # holy crap I dodged a bullet here. 
                # Silva will include numbers after a genus when the OTUs of a set of seqeunces cluster into groups, 
                # but there's no formal taxonomy that would name these smaller groups into new genera. They're noted 
                # as clustering into several sub-genera by the classifier (RDP, typically), not IJSEM. 
                # It's not official taxonomy, as in it's not a name in standing in the literature and IJSEM hasn't recognized it. 
                # It can be ignored. 
                # see Henderson G, PeerJ. 2019 Mar 5;7:e6496. 

                elif len(lineitems[1].split(' ')) > 1:
                    genus=lineitems[1].split(' ')[0]
                else:
                    genus=lineitems[1]
            if lineitems[2]=='species':
                # ncbi and silva include the genus in the species taxon
                # gg does not. GG usually doesn't have any species information. 
                # Just handle that here and get rid of the -db option
                
                # first check for a blank
                if not lineitems[1]:
                    # not real elegant, but works for now
                    check_species=["unavailable", "unavailable"]
                else:
                    check_species=lineitems[1].split(' ')
                    
                # check for any bad words
                if set(bad_words).intersection(set(check_species)):
                    species="unavailable"
                # is there only one element
                elif len(check_species)==1:
                    # if genus was set to unavailable, then the species is unavailable
                    if genus=="unavailable":
                        species="unavailable"
                    else:
                        species="{} {}".format(genus, lineitems[1])
                # if it's got this far, I think it's fine
                # but it could have extra information
                # like 'Streptococcus anginosus SK52 = DSM 20563'
                else:
                    if args.verbose:
                        print("hello from validate match - {}".format(lineitems[1]))
                    species=lineitems[1]
                    
                # then add to dictionary
                results_dict[lineitems[0]]={'taxon':species, 'confidence':float(lineitems[3])}
                        
    if args.debug:
        print(results_dict)

    validate=list()

    for k,v in translate_dataset.items():
        get_map_genus_species = v.split(' ')
        # just use the first two elements
        map_genus_species = " ".join([get_map_genus_species[0], get_map_genus_species[1]])
        
        if args.verbose:
            print("\tfrom the map file:\n\t\tget_map_genus_species={}\n\t\tmap_genus_species={}\n\tfrom the query dictionary\n\t\tk={}\n\t\tv={}\n".format(get_map_genus_species, map_genus_species,k,v))
            
        # compare the names, but first check that the result is only genus:species
        get_results_genus_species=results_dict[k]["taxon"].split(' ')
        if len(get_results_genus_species) > 2:
            results_genus_species="{} {}".format(get_results_genus_species[0], get_results_genus_species[1])
        else:
            results_genus_species=results_dict[k]["taxon"]
            
        if args.verbose:
            print("\tfrom the results dictionary:\n\t\tget_results_genus_species={}\n\t\tresults_genus_species={}\n\tfrom the query dictionary\n\t\tk={}\n\t\tv={}\n".format(get_results_genus_species, results_genus_species,k,v))
          
            
        # if the names match, add them to the validate list as a match
        if map_genus_species == results_genus_species:
            validate.append("{0},{1},{2},{3},{4}".format(k, v, results_genus_species, results_dict[k]["confidence"],1))
        # if not, there's more work
        else:
            # are there known synonyms of the query results?
            if results_genus_species in synonym_dict.keys(): 
                if results_genus_species!=synonym_dict[results_genus_species]: 
                    if args.debug:
                        print("\tmap_genus_species ({0}) != results_genus_species ({1}), and {1} != the linked current name '{2}',\n\tso checking if '{0}'== '{2}'\n".format(map_genus_species, results_genus_species, synonym_dict[results_genus_species]))
                        
                    print("\tthe result name '{0}' was found to be a synonym linked to the current name '{1}'\n\tcomparing the mapped name '{2}' to the current name '{1}'".format(results_genus_species, synonym_dict[results_genus_species], map_genus_species))
                    
                    if map_genus_species==synonym_dict[results_genus_species]:
                        print("\tand found to be the same, writing [{0},{1},{2},{3},{4}] to file\n".format(k, v, synonym_dict[results_genus_species],results_dict[k]["confidence"],1))
                        validate.append("{0},{1},{2},{3},{4}".format(k, v, synonym_dict[results_genus_species], results_dict[k]["confidence"],1))
                    else:
                        print("\tand found to be different, writing [{0},{1},{2},{3},{4}] to file\n".format(k, v, synonym_dict[results_genus_species],results_dict[k]["confidence"],0))
                        validate.append("{0},{1},{2},{3},{4}".format(k, v, synonym_dict[results_genus_species], results_dict[k]["confidence"],0))
                else:
                    if args.debug:
                        print("\tresults_genus_species ({0}) == synonym_dict[results_genus_species] ({1})\n\tand map_genus_species ({2}) != {0}, so writing to file as a mismatch\n".format(results_genus_species, synonym_dict[results_genus_species], map_genus_species))
                    validate.append("{0},{1},{2},{3},{4}".format(k, v, synonym_dict[results_genus_species], results_dict[k]["confidence"],0))
                    
                                                 
            # no match from above, and no synonyms, means there is no match 
            # however, there could be misspellings. Something for later.
            else:
                if args.debug:
                    print("\tmap_genus_species != results_genus_species and no synonyms found means there is no match. However, there could be misspellings. Something for later.\n")
                validate.append("{0},{1},{2},{3},{4}".format(k, v, results_genus_species, results_dict[k]["confidence"],0))

        
    print("writing to {}".format(outfile))
    with open(outfile, 'w') as g:
        g.write("id,query,blca,confidence,match\n")
        for x in validate:
            g.write(x+'\n')








