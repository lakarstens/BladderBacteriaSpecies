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
    blcaout_to_csv converts the outfile formatted by BLCA to a csv that can be used as a dataframe.
    
    The positional argument is the path to a folder, and the files in the folder are either different
    sequences to be processed separately, or sequences that are to be consolidated into one (as in, a 
    multirecord fasta file was split up for parallel processing)
    
    3.9.2020 update
    Queries labeled as 'unclassified' were previously filled out with rank=NA and confidence=0. 
    This is bad. Unclassified needs to be put in the missed match box (cell C), so keep the designation. 
    'Unavailable' is reserved for when there's no taxon information at a particular rank.
       
    
    Use the --consolidate option to consolidate files, otherwise the files are processed seperately
    
""", epilog="This script has no other dependencies")

parsing.add_argument("--infolder", "-i", required=True, help="The full or relative path to the directory. How the files are treated depends on the --consolidate option.")
parsing.add_argument("--consolidate", "-c", action="store_true", help="Write all files in directory as one output file")
parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
  
args=parsing.parse_args()

#--------------------------
#   get to work
#--------------------------

def process_files(dict_of_file, outfile):
    second_dict=dict()
    for k,v in dict_of_file.items():
        if args.verbose:
            print("key={} value={}".format(k,v))
        
        # split from phylum:Proteobacteria;100.0
        # into rank:taxon, confidence level
        taxons=v.split(";")  
            
        # if the sequence can be classified
        if v!="Unclassified":
            if taxons[12].split(":")[1]=="":
                species_name="unavailable"
            else:
                species_name=taxons[12].split(":")[1]
                
            second_dict[k]={"domain":(taxons[0].split(":")[1], taxons[1]), "phylum":(taxons[2].split(":")[1], taxons[3]), "class":(taxons[4].split(":")[1], taxons[5]), "order":(taxons[6].split(":")[1], taxons[7]), "family":(taxons[8].split(":")[1], taxons[9]), "genus":(taxons[10].split(":")[1], taxons[11]), "species":(species_name, taxons[13])}
            
        # if not, fill with NA and 0 <- not anymore
        # 3/9 - the update
        else:
            second_dict[k]={"domain":("unclassified", 0), "phylum":("unclassified", 0), "class":("unclassified", 0), "order":("unclassified", 0), "family":("unclassified", 0), "genus":("unclassified", 0), "species":("unclassified", 0)}

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

# write all infiles to one outfile
if args.consolidate:       
    outfile="../processed_files/formatted_{1}_{0}.csv".format(runtime, infile.split("/")[-1])
    
    if args.verbose:
        print("\tlooking for files in folder {}".format(infile))
           
    first_dict=dict()
    for gf in get_files:           
        with open("{0}/{1}".format(infile,gf), 'r') as f:
            r=csv.reader(f, delimiter='\t')
            for row in r:
                first_dict[row[0]]=row[1]
                
    process_files(first_dict, outfile)

# write each infile to seperate outfiles
else:
    print("writing to folder {}".format(outfolder))
    for gf in get_files:
        outfile="{2}/formatted_{1}_{0}.csv".format(runtime, gf.split("/")[-1].split(".")[0], outfolder)
            
        if args.verbose:
            print("\treading lines in file {}".format(gf))
                
        first_dict=dict()
        with open("{}/{}".format(infile,gf), 'r') as f:
            r=csv.reader(f, delimiter='\t')
            for row in r:
                first_dict[row[0]]=row[1]
                
        process_files(first_dict, outfile)


    """
    key=JGZX01000016.1066724.1068272 value:
            superkingdom:Bacteria
            100.0
            phylum:Proteobacteria
            100.0
            class:Gammaproteobacteria
            100.0
            order:Pseudomonadales
            100.0
            family:Pseudomonadaceae
            100.0
            genus:Pseudomonas
            88.9465534466
            species:Pseudomonas aeruginosa
            36.431862581
            
    'FNRU01000002.917551.919071': {
     'domain': ('Bacteria', '100.0'),
     'phylum': ('Actinobacteria', '100.0'),
     'class': ('Actinobacteria', '100.0'),
     'order': ('Corynebacteriales', '100.0'),
     'family': ('Corynebacteriaceae', '100.0'),
     'genus': ('Corynebacterium', '100.0'),
     'species': ('Corynebacterium pseudotuberculosis', '28.8920912421')}
     
    id              taxon                   rank    confidence
    GU812308.1.1597 Bacteria                domain  100.0
    GU812308.1.1597 Firmicutes              phylum  100.0
    GU812308.1.1597 Bacilli                 class   100.0
    GU812308.1.1597 Lactobacillales         order   100.0
    GU812308.1.1597 Lactobacillaceae        family  100.0
    GU812308.1.1597 Lactobacillus           genus   100.0
    GU812308.1.1597 Lactobacillus casei     species 39.5682983683
    """
            
#print("id\ttaxon\trank\tconfidence")



    
    
    
    