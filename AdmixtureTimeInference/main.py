"""This is the main module for AdmixtureTimeInference project

This module contains two functions: get_tract_lengths finds homozygosity tracts in simulated data and outputs tract lengths, while get_mle inference time since admixture from those tract lengths
"""

import file_handling as flh
import parsers as prs
import numpy as np 
import pandas as pd
import pickle
import argparse
import re 


def get_tract_lengths(input_file, mapfile, admixture_type, sims=[0,100],drop_ends=False):
    '''
    Takes in an input bp file and returns lists of lists ( one for each simulation in bp file ) continaing tract (homozygous tract) lengths ( in units of genetic position ). Returns a dictionary with population identifiers as keys continaing the lists of lists. 
    based on admixture type (list of founder populations). sims allows you to choose range of samples to return. drop_ends optionsallows you to specify if segment at end of chromosome is to be included
    '''
#{{{
    #Check if admixture_type is compatible with pedigree structure ( length of admixture_type list should be consistent with number of meiosis simulated, which is pulled using re.search from the input_file name ) 

    print('Parameters printed from main.py: input_file = {}, mapfiles = {}, sims = {} , admixture_type = {}, drop_ends = {}'.format(input_file,mapfile,sims,admixture_type,drop_ends))
            
    meios = int(re.search('_(\d)meios', input_file).group(1))
    assert 2**(meios + 1) == len(admixture_type) 

    #initialize tract lengths list
    tract_lengths_list_A = []
    tract_lengths_list_B = []

    #Get focal samples from input file 
    exp1 = 'g' + str(meios+2) 
    bp_list_focal = flh.print_all_match(input_file, exp1) #List containing all strings corresponding to focal sample given meios

    bp_list_subset = bp_list_focal[sims[0]*2:sims[1]*2] #get simulations corresponding to input

    #Convert physical positions to interpolated genetic positions 
    bp_list_subset_genetic = prs.phys2gen_list(bp_list_subset, mapfile) 

    #Loop over the simulations
    for i in range(len(bp_list_subset_genetic)):
        if i%2 != 0: #we want to loop in batches of two to process two haplotypes together
            continue

        #Ancestry parse the input strings corresponding to the current simulation. If simulation is weird and 
        #doesn't contain the right number of unique haplotypes, ancestry_parse returns an empty list and we continue
        focal_ancestry_bp = prs.ancestryParse([bp_list_subset_genetic[i],bp_list_subset_genetic[i+1]],admixture_type)

        #Get dictionary of homozygous tracts for the current simulation 
        if drop_ends == True:
            focal_ancestry_parsed = [prs.haplo_parse(i,drop_ends=True) for i in focal_ancestry_bp]
        else:
            focal_ancestry_parsed = [prs.haplo_parse(i) for i in focal_ancestry_bp]

        focal_ancestry_homoz_tracts = [prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'A'),\
                prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'B')] 

        #pull out tract lengths for each population, append to corresponding list, return lists 
        tract_lengths_A = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[0]] 
        tract_lengths_list_A.append(tract_lengths_A)
        tract_lengths_B = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[1]] 
        tract_lengths_list_B.append(tract_lengths_B)
    
    returnDict = {"A":tract_lengths_list_A,"B":tract_lengths_list_B}
#}}}
    return(returnDict)


def get_mle_list(input_list):
    '''
    Takes in list of lists of tract lenghts ( one for each simulation in bp file ) and 
    returns mle estimate of admixture time for each sample as a list 
    '''
#{{{
    mle_list = []

    for sample_list in input_list:
        mle_sample = len(sample_list)/(sum(sample_list)*2)
        mle_list.append(mle_sample)    
#}}} 
    return(mle_list)


if __name__ == "__main__":

    # initiate the parser
    parser = argparse.ArgumentParser()

    # add long and short argument
    parser.add_argument("--inputfile", "-i",required=True, help="set input file (must be a bp file)")
    parser.add_argument("--mapfile", "-m", required=True, help="provide path to map file for conversion of physical positions in bp file to genetic positions")
    parser.add_argument("--admixture", "-a", required=True,help="admixture type as string: e.g 'ABAB'")
    parser.add_argument("--nsims","-s", nargs='+', required=True,help="range of simulations for which to store results: e.g 0 100")
    parser.add_argument("--drop_ends","-d", help='Specify as True if last blocks are to be dropped')
    parser.add_argument("--outfile", "-o", required=True,help="e.g outfile.pkl")

    # read arguments from the command line
    args = parser.parse_args()

    #print(args.inputfile, args.generations, args.admixture[0],type(args.admixture[1]), args.nsims[0], args.nsims[1])
     
    #~~~~~~~~ Body ~~~~~~~~~~~~~~~~~ #
    #Pre-processing user inputs before calling get_tract_lengths

    args.admixture = list(args.admixture)
    args.nsims = [int(x) for x in args.nsims]
    if args.drop_ends:
        print('dropping ends of chromosomes')
        res = get_tract_lengths(args.inputfile, args.mapfile, args.admixture, args.nsims, args.drop_ends)
    else:
        print('drop_ends argument not specified, retaining chromosome ends')
        res = get_tract_lengths(args.inputfile, args.mapfile, args.admixture, args.nsims)

    #pickle the lists for later use
    with open(args.outfile,"wb") as f:
        pickle.dump(res, f)
