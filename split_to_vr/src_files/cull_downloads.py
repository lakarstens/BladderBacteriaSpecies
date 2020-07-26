#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime
from Bio import SeqIO

#--------------------------
#       argparse!
#--------------------------

parsing=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    The first time I gathered all the bladder sequences, I didn't pay real close attention to what Silva was giving me. 
    The download results are described in the file thesis/gather_markers/ssu16s/bladder_population/silva_download4_20.txt. 
    It's not wrong, but there were other accession numbers of the same species that I searched from the KTW dataset. 
    Some of the species also included the other copies of the 16S gene that are present in the bacterial genome. 
    These may be sequences obtained from environmental samples, or grown in experiments that are not representative 
    of the bladder environment. In total, there were 150 records. I've been using this set the whole time, but now 
    I realize that I should be using the type specimens that are deposited into ATCC or other holding facilities when 
    the bacterial species is given a name in good standing. That way, the sequences that I'm using are the benchmark 
    to what environmental samples can be compared.
    """, epilog="no other dependencies except biopython")

parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
parsing.add_argument("--infile", "-i", required=True, help="the path to the multi-record fasta file downloaded from Silva")
  
args=parsing.parse_args()

#--------------------------
#   get to work
#--------------------------


runtime=datetime.now().strftime('%Y-%m-%d_%M')

with open("../../resources/accessions.txt", 'r') as f:
    acc_info=f.readlines()
    
accessions=[x.split('|')[4].strip() for x in acc_info]

#get_records={x.id.split('.')[0]:x for x in SeqIO.parse("arb-silva.de_2019-11-10_id739827_tax_silva.fasta", "fasta")}
get_records={x.id.split('.')[0]:x for x in SeqIO.parse(args.infile, "fasta")}

get_ids=[l for l in get_records.keys()]

smaller=[g for g in get_ids if g in accessions]

#print(len(smaller))

print("writing to ktw_16s_type_{}.fna".format(runtime))
with open("ktw_16s_type_{}.fna".format(runtime), 'a') as o:
    for s in smaller:
        SeqIO.write(get_records[s], o, 'fasta')
    