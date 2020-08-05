#!/usr/bin/env python3

from Bio import AlignIO
import sys
import swa_tools
from datetime import datetime
import os
import argparse


# get file path and load MSA
# make sure there's a file included

#--------------------------
#       argparse!
#--------------------------

parsing=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    sized_window() takes the text file "infile" as a positional argument. 
    Each line in "infile" holds a path to a multisequence alignment and a list of optimal window sizes separated by commas. Like this: 

        ../raw_data/bladder_16s_output/bladder_16s.aln 20,70,140,170
        
    The script is run like this:
           
        python sized_window.py /path/to/infile.txt

    The output is written to a new folder with the format:

            ../processed_data/specific_windows_weighted_<runtime>
        
""", epilog="This script depends on the following swa_tools functions:\n  sliding_window()\n{}".format(swa_tools.sliding_window.__doc__))
    
parsing.add_argument("infile", help="The full or relative path to a batch file.")    
parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
parsing.add_argument("--debug", "-d", action="store_true", help="more detail about what is happening")

args=parsing.parse_args()

#--------------------------
#       get to work
#--------------------------

runtime=datetime.now().strftime('%m_%d_%H%M')
new_dir="../processed_data/specific_windows_weighted_{0}".format(runtime)
if not os.path.exists(new_dir):
    os.makedirs(new_dir)
    if args.verbose:
        print("\tcreating new directory in {}".format(new_dir))
   

#print("\twriting to folder {}".format(new_dir))

if args.verbose:
    print("\tlooping through lines of infile")
    
#---------------------------------------
#   create the weighted entropy file
#---------------------------------------

#alignment_files=os.listdir(args.infile)
#alignment_files_path = line.split()

alignment_path=list()

with open(args.infile, 'r') as h:
    for line in h:
        alignment_lines=line.split()
        alignment_path.append(alignment_lines[0])
        
if args.debug:
    print(alignment_path)

for af in alignment_path:
    # read the MSA clustal file format with AlignIO
    alignment = AlignIO.read("{}".format(af),"clustal")
    
    if args.debug:
        print("\tsending {} to AlignIO".format(af))

    if args.verbose:
        print("\tworking on {}".format(af))
        
    weighted_dict=swa_tools.create_weight_dict(alignment)
    
    #if args.verbose or args.debug:
     #   print("\n\tsending {} dictionary to write_weighted_ent()".format(weighted_dict))
     
    get_base = os.path.basename(af)
    basename = os.path.splitext(get_base)[0]
    
    if args.debug:
        print("\n\tsending\n\t\t{}\n\t\t{}\n\t\t{}\n\tto write_weighted_ent() function".format(alignment, weighted_dict, basename))
        
    path_to_weight_values = swa_tools.write_weighted_ent(alignment, weighted_dict, basename)
    
    if args.debug:
        print(path_to_weight_values)
    
#-------------------------------------------
#   loop through different sized windows
#-------------------------------------------

window_dict={}  

if args.verbose or args.debug:
    print("\tcreating dictionary of source file and window sizes from {}".format(args.infile))

 
with open(args.infile, 'r') as f:
    for line in f:
        cog_windows=line.split()
        
        if args.debug:
            print(cog_windows)
            
        window_dict[cog_windows[0]] = [int(x) for x in cog_windows[1].split(',')]

#track_progress=0
if args.verbose or args.debug:
    print("\twindow_dict = {}".format(window_dict))
    
for k,v in window_dict.items():
    # probably better to make the file argument contain full paths 
    #alignment = AlignIO.read("../raw_data/aln/{0}".format(k),"clustal")
    #alignment = AlignIO.read(k,"clustal")
    get_filename=k.split("/")
    file_prefix=get_filename[-1][0:7]
    if args.verbose:
        print("\tworking on\n\t\tk={}\n\t\tv={}".format(k, v))
        
    #progress=track_progress/float(len(k))
    #swa_tools.update_progress(progress, "Sliding windows across {0}... ".format(file_prefix)) 

    for n in v:   
        if args.verbose or args.debug:
            print("\tsending \n\t\tn={}\n\t\tk={}\n\t\tfile prefix={}\n\t\tnew dir={} to sliding_window()".format(n, k, file_prefix, new_dir))
            
        # changed k to the calculated weighted entropy scores for the msa
        swa_tools.sliding_window(n, path_to_weight_values, file_prefix, new_dir)
        
    #track_progress+=1
   
    
    



    
    

