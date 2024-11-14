#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:22:25 2024

@author: avierstr
"""

import edlib
import re
import argparse
import os
import sys
from Bio import SeqIO
import csv
import copy
#==============================================================================
version = '2024-11-13'  # version of the script
#==============================================================================
def get_arguments():
    def range_limited_float_type(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 0 or f > 1:
            raise argparse.ArgumentTypeError("Argument must be > " + str(0) 
                                             + " and < " + str(1))
        return f
    def valid_file(param): 
        paramlist = [] # make list of files to process
        if os.path.isfile(param): # if input is file
            # check if input file ends on .fastq or .fasta
            base, ext = os.path.splitext(param)
            if ext.lower() not in ('.fasta', '.fastq'): 
                raise argparse.ArgumentTypeError('File extension must be .fastq or .fasta') 
            paramlist.append(param)
        elif os.path.isdir(param): # if input is folder 
            with os.scandir(param) as iterator:
                for file in iterator:
                    if file.name.lower().endswith('.fastq') or file.name.lower().endswith('.fasta'):
                        paramlist.append(file.path)
            paramlist.sort()
            if len(paramlist) == 0:
                sys.exit('Can not find files in folder to process.  File extension must be .fastq or .fasta')
        else:
            sys.exit('Can not find a file or folder to process.  File extension must be .fastq or .fasta')
        param = paramlist
        return param
    
    def dir_path(string):
        string = os.path.join(os.getcwd(), string)
        if not os.path.exists(string):
            os.makedirs(string) # create the folder
        return string

    parser = argparse.ArgumentParser(description='map_my_primers, a tool to map primers on sequences' )
    parser.add_argument('-i', '--input', required=True, default = os.getcwd(), type = valid_file,
                        help='Input folder or file in fastq or fasta format')
    parser.add_argument('-o', '--outputfolder', type=dir_path, 
                         help='Save the results in the specified\
                            outputfolder. Default = mapped')
    parser.add_argument('-min', '--minlength', type = int, default=1,
                        help='Minimum readlenght to process.')
    parser.add_argument('-max', '--maxlength', type = int, 
                        help='Maximum readlenght to process.  Default = No limit')
    parser.add_argument('-pr', '--primers', type = str, required=True,
                        help='File with primers/barcodes used in experiment.')
    parser.add_argument('-er', '--error', type = range_limited_float_type,
                        default=0.15, help='Percentage error allowed in editdistance for\
                            adapters and barcodes. Default = 0.15 (15 percent)')
    args = parser.parse_args()
    # print(args)
    # sys.exit()
    return args
#==============================================================================
def align(q, t, m, a, k):
    s = edlib.align(q, t, mode=m, task=a, k=k, # task='path' locations
                    additionalEqualities=[("R", "A"), ("R", "G"),
                          ("Y", "C"), ("Y", "T"),
                          ("M", "A"), ("M", "C"),
                          ("K", "G"), ("K", "T"),
                          ("S", "G"), ("S", "C"),
                          ("W", "A"), ("W", "T"),
                          ("N", "A"), ("N", "T"), 
                          ("N", "G"), ("N", "C")])
    return s
#==============================================================================
def compl_reverse(self):
    '''
    make the complement reverse of the sequence
    '''
    inp  = 'ATCGRYKMSWN' # translate table for complement
    outp = 'TAGCYRMKSWN'
    complement = ''.maketrans(inp, outp)
    R = (self[::-1]).translate(complement)  # complement reverse
    return R
#============================================================================== 
def put_in_dict(key, value, dic):
    if key in dic:
        adapt = dic.get(key)
        adapt.append(value)
        dic[key] = adapt
    else:
        dic[key] = [value]
    return dic
#============================================================================== 
def put_in_dict2(key, value, dic):
    if key in dic:
        adapt = dic.get(key)
        adapt.append(value)
        dic[key] = adapt
    else:
        if type(value) is list:
            dic[key] = value
        else:
            dic[key] = [value]
    return dic
#==============================================================================  
def read_custom_barcodes(file):
    ''' 
    Read the file with custom barcodes and put them in list.  If "forward" and "reverse" 
    is used in the name, it means they are different, otherwise they are the same.
    '''
    try:
        used_bc_list = []
        # with open(os.path.join(infolder, file), newline='') as csvfile:
        with open(file, newline='') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read()) # check if comma or tab separated
            csvfile.seek(0) # got back to begin of file
            reader = csv.reader(csvfile, dialect)
            for row in reader:
                if any(x.strip() for x in row): # remove empty lines
                    if not row[0].startswith('#'): # comment in line
                        name = row[0]
                        fbc = row[1].strip().replace(' ','').upper()
                        used_bc_list.append([name, fbc])
                        name = row[0]
                        rbc = row[2].strip().replace(' ', '').upper()
                        used_bc_list.append([name, rbc])
    except FileNotFoundError:
        print('Can not find primers file ' + file)
        sys.exit()
    return used_bc_list
#==============================================================================
def create_search_list(used_bc_list):
    '''
    All bc (also merged ones) are 5'-3' direction. 
    front_bc is for barcodes on the 5' side of each read
    rear_bc is for the complement reverse barcodes on the 3' side of the read.
    Nested lists: [[name, F_bc, R_bc],...]  [[name, F_cr_bc, R_cr_bc],...]
    '''
    front_bc = sorted(used_bc_list, key=lambda x: x[0])
    front_bc = [x for x in front_bc if len(x[1]) > 0]
    rear_bc = [[name, compl_reverse(seq)] for name, seq in front_bc]
    return front_bc, rear_bc
#==============================================================================
def filter_locations(locations):
    """
    check if there are multiple locations where the adapter fits at the same
    start position and only use the first for each start position.
    """
    if len(locations) > 1:  
        for l in range(0, len(locations)-1):
            for m in range(l+1, len(locations)):
                if locations[l][0] == locations[m][0]:
                    locations[m] = ('','')
        locations = [x for x in locations if x != ('','')]
    return locations
#============================================================================== 
def best_hits(hitlist):
    """
    filter for the best and longest hits for each position
    """
    hitlist2 = []
    # if there are several locations for a hit, duplicate them per location
    # ['16S', 'AAGTCGTAACAAGGTAACCGTCATCGACCATC', [('editDistance', 0), ('alphabetLength', 4), ('locations', [(1519, 1550), (3076, 3107)]), ('cigar', '32=')], 1]
    for x in hitlist:
        loc = x[2][2][1]
        # hitlist.pop(i)
        for j in range(len(loc)):
            x[2][2] = list(x[2][2]) # tuple can not be changed
            x[2][2][1] = [loc[j]]
            x[2][2] = tuple(x[2][2])
            hitlist2.append(copy.deepcopy(x))
    hitlist = hitlist2    
    
    hitlist.sort(key=lambda x: x[2][2][1][0][1]) # sort based on location

    for x in range(0, len(hitlist)-1):
        for y in range(x+1, len(hitlist)):
            dis1, loc1 = hitlist[x][2][0][1], hitlist[x][2][2][1][0]
            dis2, loc2 = hitlist[y][2][0][1], hitlist[y][2][2][1][0]
            loc = set(loc1)
            loc.update(loc2)
            if len(loc) == 3: # start or stop is the same
                if dis1 == dis2: # distance is the same
                    if abs(loc1[0]-loc1[1]) > abs(loc2[0]-loc2[1]):
                        hitlist[y][-1] = ('') # take longest hit
                    else:
                        hitlist[x][-1] = ('')
                if dis1 != dis2: # distance is the different
                    len1 = abs(loc1[0]-loc1[1]) # length of the hit
                    len2 = abs(loc2[0]-loc2[1])
                    if len1 > len2:
                        if dis1/len1 < 0.15:
                            hitlist[y][-1] = ('') # take longest hit
                        else:
                            hitlist[x][-1] = ('') # take longest hit
                    elif len2 > len1:
                        if dis2/len2 < 0.15:
                            hitlist[x][-1] = ('') # take longest hit
                        else:
                            hitlist[y][-1] = ('') # take longest hit
            elif len(loc) == 2: # start and stop is the same
                if dis1 < dis2:
                    hitlist[y][-1] = ('') # remove worst hit
                elif dis1 > dis2:
                    hitlist[x][-1] = ('') # remove worst hit
                elif dis1 == dis2:
                    hitlist[y][-1] = ('') # remove one hit
            elif len(loc) == 4: # start and stop are different
                a = range(loc1[0], loc1[1])
                b = range(loc2[0], loc2[1])
                if set(a).issubset(b): # if a is part of b
                    hitlist[x][-1] = ('') # remove shortest hit
                elif set(b).issubset(a): # if b is part of a
                    hitlist[y][-1] = ('') # remove shortest hit

    hitlist = [x for x in hitlist if x[-1] != '']
    
    hitlist = [str(x) for x in hitlist] # make string of list items
    hitlist = list(set(hitlist))
    hitlist = [eval(x) for x in hitlist] # make list of string items
    hitlist.sort(key=lambda x: x[2][2][1][0][1]) # sort based on location

    return hitlist
#==============================================================================
def find_hits(rec, front_bc, rear_bc, error):
    number = re.compile(r'\d+') # the numbers to find in edlib result
    symbol = re.compile(r'\D') # the letters or symbols to find in edlib result
    m = 'HW'
    a = 'path'
    hitlist = [' ']
    hitlist2 = []
    search_list = [x for x in [front_bc, rear_bc] if len(x) > 0]
    # for lis in search_list:
    #     lis.sort(key=lambda x: len(x[1]), reverse=True) 
    seq = [x for x in rec.seq] # make a list of the sequence
    seq2 = seq[:]
    alignlist = [' ']*len(seq)
    namelist = [' ']*len(seq)
    """
    if a 100% hit is found and there is also a 99% hit, it will only show the 100% hit
    for that reason the script has to scan 2 times
    """
    k1 = -1 # minimum error
    while len(hitlist) > 0:
        hitlist = []
        for i, lis in enumerate(search_list):
            for name, BC in lis:
                BC = [x for x in BC]
                k = len(BC)*error
                s = align(BC, seq2, m, a, k) 
                if k > (s['editDistance']) > k1: # if a approximate hit is found
                    # print(s['locations'])
                    loc = filter_locations(s['locations'])
                    s['locations'] = loc
                    # print(s)
                    hitlist.append([name, ''.join(BC), list(s.items()), i])
        
        hitlist = best_hits(hitlist)
        
        for h in hitlist:
            hitlist2.append(h)
            BC = [x for x in h[1]]
            k = len(BC)*error
            s = align(BC, seq2, m, a, k) 
            # recalculate alignments because gaps can influence locations
            if k > (s['editDistance']) > k1: # if a approximate hit is found
                loc = filter_locations(s['locations'])
                for j in range(len(loc)):
                    begin, end = s['locations'][j]
                    seq2[begin:end] = '-'*(end-begin) # remove the hit zone from the intermediate sequence
                    nu = re.findall(number, s['cigar'])
                    sy = re.findall(symbol, s['cigar']) 
                    scorelist = []
                    for x, y in zip(nu, sy):
                        scorelist += int(x) * [y]
                    for i, x in enumerate(scorelist):
                        if x =='D':
                            BC.insert(i, '-') # insert gap in BC
                        if x == 'I':
                            seq.insert(begin + i, '-') # inser gap in seq   
                            seq2.insert(begin + i, '-') # inser gap in seq   
                            namelist.insert(begin + i, ' ') # inser gap in namelist   
                    alignlist[begin:end+1] = BC
                    if h[-1] == 0:
                        name = 'FW_' + h[0] + '>'
                        endname = len(name)
                        namelist[(end-endname)+1:end+1] = name
                    elif h[-1] == 1:
                        name = '<RV_' + h[0]
                        endname = len(name)
                        namelist[begin:endname+1] = name
        k1 += 1  
    with open(outfile, 'a') as of:
        print('>' + str(rec.description), file=of)
        for x in hitlist2:
            print(x[0],x[2], file=of)
        for i in range(0, len(seq), 100):
            print(''.join(seq[i:i+100]), file=of)
            print(''.join(alignlist[i:i+100]), file=of)
            print(''.join(namelist[i:i+100]), file=of)
#==============================================================================
if __name__ == '__main__':
    args = get_arguments()
    outputfolder = args.outputfolder
    infolder_file_list = args.input
    min_length = args.minlength
    max_length = args.maxlength
    primers = args.primers
    error = args.error
    used_bc = read_custom_barcodes(primers)
    front_bc, rear_bc = create_search_list(used_bc)
    
    for infolder_file in infolder_file_list:
        infolder, infile = os.path.split(os.path.realpath(infolder_file))
        # outfolder, ext = os.path.splitext(infile)
        if not outputfolder: # if outputfolder is not given
            try:
                outputfolder = os.path.join(infolder, 'mapped')
                os.makedirs(outputfolder)
            except FileExistsError:
                pass

        with open(os.path.join(infolder, infile), 'r') as inf: # check the fileformat
            line = inf.readline()
            if line[0] == '>':
                fileformat = 'fasta'
            elif line[0] == '@':
                fileformat = 'fastq'
        
        outfile = os.path.join(outputfolder, infile.replace('.fasta', '_map_primers.fasta').
                               replace('.fastq', '_map_primers.fasta'))
        try: # remove outfile if exists
            os.remove(outfile)
        except FileNotFoundError:
            pass
        
        q = 0
        with open(os.path.join(infolder, infile), "r") as handle: 
            print(infile)
            for record in SeqIO.parse(handle, fileformat):
                record = record.upper() # make all sequences uppercase
                if len(record.seq) >= min_length:
                    if max_length is None:
                        rec = record
                        find_hits(rec, front_bc, rear_bc, error)
                        q += 1
                    else: 
                        if len(record.seq) <= max_length:
                            rec = record
                            q += 1
                            find_hits(rec, front_bc, rear_bc, error)
                print('\rNumber of reads processed: ' + str(q), end='')
        print('')