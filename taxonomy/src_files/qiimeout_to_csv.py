#!/usr/bin/env python3

import csv
import os
from datetime import datetime
import sys
import argparse

#infile="test_16s_v4_blcaout.txt"
#infile="../raw_data/test_16s_v4.fna.blca.out"

#--------------------------
#       argparse!
#--------------------------

parsing=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    qiimeout_to_csv converts the outfile formatted by Qiime to a csv that can be used as a dataframe.
    
    if genus or species taxons from the qiime outfile are blank e.g. 'g__; s__' or some combination, then the taxon is labeled 'unavailable'. 
    I thought if both are blank I should call them 'unclassified', but I'm going to leave it up to the classifer 
    to make up it's own mind about whether or not it can classify a query.
    
    The positional argument is the path to a folder, and the files in the folder are different sequences to be processed separately
    
    5/19 update
    Also used a pretrained database, plenty of problems with formatting to capture.
        e.g. the ranks are named "D_5__" instead of "g__". Added a clause to catch that.
        no space after the ";", had to catch that
        etc

    
""", epilog="This script has no other dependencies besides BioPython")

parsing.add_argument("--infolder", "-i", required=True, help="The full or relative path to the directory. How the files are treated depends on the --consolidate option.")
parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
parsing.add_argument("--debug", action="store_true", help="print the nitty-gritty for debugging purposes")
  
args=parsing.parse_args()

#--------------------------
#   get to work
#--------------------------

def process_files(dict_of_file, outfile):
    second_dict=dict()
    for k,v in dict_of_file.items():
        if args.verbose:
            print("key={} value={}".format(k,v))
        
        # split from 'Y17005.1.1451': {'tax_chain': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Aerococcaceae; g__Aerococcus; s__christensenii', 'confidence': 0.9999094074904588}}
        # into rank:taxon, confidence level
        get_taxchain=[x.strip() for x in v['tax_chain'].split(";")]  
        
        if args.verbose:
            print(get_taxchain)
        
        taxons=list()
        # need to put in a clause that adds 'unavailable' for blank taxons
        for x in get_taxchain:
            split_chain=x.split('__')
            if (split_chain[0]=='g' or split_chain[0]=='D_5'):
                if split_chain[1]=="":
                    genus='unavailable'
                else:
                    genus=split_chain[1]   
                taxons.append((genus, 'NA'))     
                
            elif (split_chain[0]=='s' or split_chain[0]=='D_6'):
                if split_chain[1]=="":
                    species="{} {}".format(genus, 'unavailable')
                else:     
                    check_for_genus = split_chain[1].split()
                    if check_for_genus[0] == genus:
                        species="{} {}".format(genus, check_for_genus[1])
                    else:
                        species="{} {}".format(genus, split_chain[1])
                taxons.append((species, v['confidence']))
                
            else:
                taxons.append((split_chain[1], 'NA'))
                
        if args.verbose:
            print(taxons)
            
        second_dict[k]={"domain":(taxons[0][0], taxons[0][1]), "phylum":(taxons[1][0], taxons[1][1]), "class":(taxons[2][0], taxons[2][1]), "order":(taxons[3][0], taxons[3][1]), "family":(taxons[4][0], taxons[4][1]), "genus":(taxons[5][0], taxons[5][1]), "species":(taxons[6][0], taxons[6][1])}

    # write to a new csv file     
    print("writing to {}".format(outfile))
    with open(outfile, 'a') as f:
        f.write("id,taxon,rank,confidence\n")
        for k,v in second_dict.items():
            for m,n in v.items():
                f.write("{},{},{},{}\n".format(k, n[0], m, n[1]))
    
infile=args.infolder
#print(infile,infile.split("/")[-1].split(".")[0])
runtime=datetime.now().strftime('%Y-%m-%d_%H%M')

get_files=os.listdir(infile)

outfolder="../processed_files/formatted_outfiles_{}".format(runtime)
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

if args.verbose:
    print(", ".join(get_files))


print("writing to folder {}".format(outfolder))
for gf in get_files:
    outfile="{2}/formatted_{1}_{0}.csv".format(runtime, gf.split("/")[-1].split(".")[0], outfolder)
        
    if (args.verbose or args.debug):
        print("\treading lines in file {}".format(gf))
            
    first_dict=dict()
    with open("{}/{}".format(infile,gf), 'r') as f:
        r=csv.reader(f, delimiter='\t')
        for row in r:
            if row!=['Feature ID', 'Taxon', 'Confidence']:
                first_dict[row[0]]={'tax_chain':row[1], 'confidence': round(float(row[2]), 5)}
    
    if args.debug:
        print(first_dict)
        
    process_files(first_dict, outfile)


    
    
    
    