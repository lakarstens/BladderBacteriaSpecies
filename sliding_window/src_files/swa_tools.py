#!/usr/bin/env python3

from Bio import AlignIO
import math
import statistics
import sys
from datetime import datetime
import os


# read the MSA clustal file format with AlignIO
#alignment = AlignIO.read("../raw_data/tcoffee_msa_16s.aln","clustal")
#alignment = AlignIO.read(sys.argv[1],"clustal")

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

def shannon_entropy(col, verbose=False):
    """
        calculates and returns the shannon entropy of the col argument.
    
        there's probably a way to do a block of columns, 
        but AlignIO seems to insist on including the 
        record id in the multi-column slice
    """

    alpha=list(set(col))
    shannon=list()
    for a in alpha:
        freq = col.count(a)/len(col)
        # math.log() is base e by default
        logf = math.log(freq)
        shannon.append(freq*logf)
        if (verbose):
            print("\tnumber of {0} is {1}\n\t\tfrequency = {2}\n\t\tlog of frequency = {3}".format(a, col.count(a), freq, logf))
    
    entropy=abs(sum(shannon))
    if (verbose):
        print("shannon entropy of col {} is {}\n".format(col, entropy))
    return entropy
    
def sliding_window(window, file_path, gene, in_folder, verbose=False):
    """
        rewrite of sliding window that reads the file of precalculated 
        weighted entropy scores from folder: 
        
            thesis_git/sliding_window/processed_data/weighted_ent_<date>/<gene>_weight_ent.csv
    """

    runtime=datetime.now().strftime('%m_%d_%H')
    #path_to_file="../processed_data/weighted_ent_{0}/{1}_weight_ent.csv".format(date, gene)
    
    if verbose==True:
        print("working on window size {0} for {1} ...".format(window, gene))
    
    with open(file_path, 'r') as c:
        all_values=[float(x.strip().split('\t')[1]) for x in c]

    aln_length=len(all_values)
    steps=aln_length-window

    print("\twriting window size {1} to {0}".format(in_folder, window))
    for p in range(0, steps,1): 
        holder=statistics.mean(all_values[p:p+window])
        with open("../processed_data/{0}/{2}_w{1}.csv".format(in_folder, window, gene), 'a') as d:
            d.write("{}\t{}\n".format(p, holder))
            
        if verbose:
            print("p={} window={} values={} mean={}".format(p,window,all_values[p:window],holder))
            
        
def find_vr_edges(genename, se_values, start_pos, cutoff, verbose=False):
    """
    Does these steps:
        - takes coordinates of peak and window size used to make the graph
        - slides a 12 bp window to the left and right of peak coordinate
        - stops when:
            - the sum shannon entropy of the window less than or equal to .322
            - or the ends of the sequence are reached  
    """
    
    runtime=datetime.now().strftime('%m_%d_%H_%M')
        
    window=12
    aln_length=len(se_values)
    r_steps=(aln_length-start_pos)-window
    l_steps=(start_pos)-window

    max_r_width=(start_pos+125-window)
    max_l_width=(start_pos-125-window)
    
    # the right side
    r_greeting="\n{1}, peak {0}\n------------------------\n  searching for the right edge".format(start_pos, genename)
    if verbose==True:
        print(r_greeting)
            
    pos=start_pos
    for p in range(1, r_steps):     
        holder=statistics.mean(se_values[pos:pos+window])

        if ((holder<=cutoff) or (pos > max_r_width)):
            r_stopped_work="\n\tstopped here.\n\t\tstarting from {3}\n\t\tend of seq is {4}\n\t\tright edge at {2}\n\t\tstopped at position = {0}\n\t\tentropy = {1}".format(pos, holder, max_r_width, start_pos, aln_length)
            if verbose==True:
                print(r_stopped_work)
                
            break

        pos+=1
         
    if holder<=cutoff:
        r_edge_state="\n\tright edge of {1} VR is {0}\n".format(pos+window, start_pos)
    else:
        r_edge_state="\n\t*** didn't find the right edge by {0}".format(pos+window)

        if pos+window==aln_length:
            r_edge_state="{0} because the end of the sequence was reached\n".format(r_edge_state)
        elif pos==max_r_width+1:
            r_edge_state="{0} because the VR would be too large\n".format(r_edge_state)
            
    if verbose==True:    
        print(r_edge_state)
    
    if holder<=cutoff:
        right_side=pos+window
    else:
        right_side=False
        
    # write all that commentary into a logfile
    with open("../processed_files/vr_width_runlog_{0}.txt".format(runtime), 'a') as c:
        c.write("{0}{1}{2}".format(r_greeting,r_stopped_work,r_edge_state))
        
    # the left side
    l_greeting="\n  searching for the left edge".format(start_pos, genename)
    if verbose==True:
        print(l_greeting)
        
    pos=start_pos
    for p in range(1, l_steps):     
        holder=statistics.mean(se_values[pos:pos+window])
        
        if ((holder<=cutoff) or (pos < max_l_width)):
            l_stopped_work="\n\tstopped here.\n\t\tstarting from {3}\n\t\tend of seq is {4}\n\t\tleft edge at {2}\n\t\tstopped at position = {0}\n\t\tentropy = {1}".format(pos, holder, max_l_width, start_pos, aln_length)
            if verbose==True:
                print(l_stopped_work)
                
            break
        else:
            l_stopped_work=""
            
        pos-=1
        
    if holder<=cutoff:
        l_edge_state="\n\tleft edge of {1} VR is {0}\n".format(pos-window, start_pos)
    else:
        l_edge_state="\n\t*** didn't find the left edge by {0}".format(pos-window)

        if pos-window==1:
            l_edge_state="{0} because the beginning of the sequence was reached\n".format(l_edge_state)
        elif pos==max_l_width-1:
            l_edge_state="{0} because the VR would be too large\n".format(l_edge_state)
            
    if verbose==True:    
        print(l_edge_state)
            
    if holder<=cutoff:
        left_side=pos-window
    else:
        left_side=False
        
    # write all that commentary into a logfile
    with open("../processed_files/vr_width_runlog_{0}.txt".format(runtime), 'a') as c:
        c.write("{0}{1}{2}".format(l_greeting,l_stopped_work,l_edge_state))
        
    if (left_side and right_side):
        #print("[{0}, {1}]".format(left_side, right_side))
        return (left_side, right_side)
    else:
        print("{2}: no variable region demarcated for peak at {0}, with threshold of {1}".format(start_pos, cutoff, genename))
        

        
def optimize_window(msa, filename, verbose=False, keep_count=True):
    """
        requires a path to a MSA file in clustal format 
        to be passed at the command line
        
        takes a window size and multisequence alignment, 
        calculates the sum of the entropy for each column 
        in the window, 
    """
    
    aln_length=msa.get_alignment_length()
    average=mean_entropy(msa, filename=filename)
    if verbose:
        print("mean of {0} MSA is {1}".format(filename, average))
    runtime=datetime.now().strftime('%m_%d_%H_%M')
    
    track_progress=10
    for window in range(10,210,10):
        opter=list()
        #print("working on window size {} ...".format(window))
        
        steps=aln_length-window
        pos=0
                
        for p in range(1, steps):     
            
            # lying Dog from stackoverflow, with the generator inside list comp
            holder=[shannon_entropy(msa[:,x+1]) for x in range(pos, window+pos)]
            mean_holder=statistics.mean(holder)
            if mean_holder>average:
                opter.append(mean_holder)
            pos+=1
            
        if keep_count==True:
            progress=track_progress/float(200)
            update_progress(progress, "Windows are sliding across "+filename)
            track_progress+=10
            
        with open("../processed_data/{0}_optimal_{1}.csv".format(runtime, filename), 'a') as c:
            c.write("{0}\t{1}\n".format(window, len(opter)))
            
            
        #print("{}\t{}".format(window, len(opter)))
    
def mean_entropy(msa, filename="testing", verbose=False):
    """
    requires an AlignIO object read from a clustal file  
    to be passed at the command line
    
    calculates and returns the mean Shannon entropy from the full 
    length of a multi-sequence alignment.
    
    Also writes to a file "<runtime>_<foo>_full_entropy.csv" where 
    <foo> is the alignment file name and <runtime> is the date and time
    """
    runtime=datetime.now().strftime('%m_%d_%H_%M')
    aln_length=msa.get_alignment_length()
    holder=[shannon_entropy(msa[:,x+1]) for x in range(1, aln_length-1)]
    
    if (verbose):
        print("working on mean\n")
        print("using MSA \n{}".format(msa))
        
    with open("../processed_data/{0}_{1}_full_entropy.csv".format(runtime, filename), 'a') as c:
        for n in range(1,aln_length):
            c.write("{}\t{}\n".format(n, shannon_entropy(msa[:,n])))
            
    if (verbose):
        print("\nmean for the full sequence \"{0}\" is {1}".format(filename, statistics.mean(holder)))
        
    return statistics.mean(holder)
    
def create_weight_dict(msa, verbose=False):
    """
        inputs
            a multisequence alignment object, as in the 
            result of reading a MSA file with AlignIO
        
        output
            a dictionary {sequence.id : weight}, which 
            maps the sequence weight to the sequence ID
    
        calculates the weight of the sequence in a MSA as 
        described in Valdar 2002 and Heinikoff & Heinikoff 1994
    """
    aln_weights={}  
    track_progress=0
    
    # number of species considered in the msa, 78 for the 16S gene
    species_num=78
    
    for s in range(0, len(msa)):   
        pos_weight=[]
        for p in range(0, msa.get_alignment_length()):
            kx=len(list(set(msa[:,p])))
            pos=msa[s].seq[p]
            nxi=msa[:,p].count(pos)
            sum_freq=1/(kx*nxi)
            pos_weight.append(sum_freq)
            if verbose:
                print("\tpos={1} kx={0} nxi={2} sum_freq={3}".format(kx,pos,nxi,sum_freq))
                
        aln_weights[msa[s].id]=1/msa.get_alignment_length()*sum(pos_weight)
        if verbose:   
            print("\ttotal weight={}".format(1/msa.get_alignment_length()*sum(pos_weight)))

        progress=track_progress/float(species_num)
        update_progress(progress, "\tCreating weighted dictionary... ")
        track_progress+=1
        
    return(aln_weights)
  
  
def write_weighted_ent(msa, w_dict, gene, conserve=False, verbose=False):
    """
        inputs
            1) a multisequence alignment
            2) the dictionary output by create_weight_dict()
            3) the gene name
        
        output
            writes the weighted entropy of one or more genes as 
            seperate files into the indicated directory 
    """
    
    runtime=datetime.now().strftime('%m_%d_%H%M')
    new_dir="../processed_data/weighted_ent_{0}".format(runtime)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    
    # the ln is supposed to be the smaller of either 
    # the alphabet or the number of sequences.
    # with 78 sequences, the alphabet is the smaller
    lamb_t=1/math.log(5)
    #print(lamb_t)
    ct=[]
    outfile="{1}/{0}_weight_ent.csv".format(gene, new_dir)

    print("\n\twriting weighted entropy values to {}".format(outfile))
    for a in range(0,msa.get_alignment_length()):
        pa=[]
        pos=msa[:,a]
        symbols=list(set(pos))
        #print("column chars={} symbol set={}".format(pos,symbols))
        gaps=0
        for y in symbols:
            gotcha=[w_dict[msa[i].id] for i,j in enumerate(pos) if j==y]
            if y=="-":
                gaps=len(gotcha)/len(pos)
                #print(gaps)
            pa.append(sum(gotcha))
            #print(sum(gotcha))
        #print(pa)
        
        shent=[m*math.log(m) for m in pa]
        w_shent=lamb_t*sum(shent)
        tot_w_ent=(1+w_shent)*(1-gaps)
        if conserve:
            ct=tot_w_ent
        else:
            ct=1-tot_w_ent
        if verbose:
            print("t={0} gaps={1} total conservation={2}".format(w_shent, gaps, tot_w_ent))
        
        
        with open(outfile, 'a') as f:
            f.write("{0}\t{1}\n".format(a, round(ct, 3)))
        
    # need to return the path of the new outfile
    return(outfile)

if __name__== "__main__":

    print("\n")
    print("+----------------------------------------------+")
    print("|                                              |")
    print("|                 swa_tools                    |")
    print("|                                              |")
    print("+----------------------------------------------+")    
    print("\nHello from swa_tools!\n\nThe following functions are available in this file\n") 
    print("* update_progress(progress, what_doing)")
    print(update_progress.__doc__)
    print("* shannon_entropy(col, verbose=False)")
    print(shannon_entropy.__doc__)
    print("* sliding_window(window, file_path, gene, in_folder, verbose=False)")
    print(sliding_window.__doc__)
    print("* find_vr_edges(genename, se_values, start_pos, cutoff, verbose=False)")
    print(find_vr_edges.__doc__)
    print("* optimize_window(msa, filename, verbose=False, keep_count=True)")
    print(optimize_window.__doc__)
    print("* mean_entropy(msa, filename='testing', verbose=False)")
    print(mean_entropy.__doc__)
    print("* create_weight_dict(msa, verbose=False)")
    print(create_weight_dict.__doc__)
    print("* write_weighted_ent(msa, w_dict, gene, conserve=False, verbose=False)")
    print(write_weighted_ent.__doc__)



    

