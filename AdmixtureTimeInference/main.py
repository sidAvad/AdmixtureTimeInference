import file_handling as flh
import parsers as prs
import numpy as np 
import pandas as pd
import pickle
import argparse

def get_tract_lengths(input_file, gens, admixture_type, sims=[0,1000], outfile=None):
    '''
    Takes in an input bp file and returns a list of lists ( one for each simulation in bp file ) continaing tract (homozygous tracts) lengths  based on number of generations and admixture type. sims allows you to choose range of samples to return
    Optionally provide outfile (.pickle) to pickle the results
    '''

    #Check if admixture_type is compatible with pedigree structure ( lenght of admixture_type list should be equal to number of founders)
    assert 2**(gens) == len(admixture_type) 

    #initialize tract lengths list
    tract_lengths_list_A = []
    tract_lengths_list_B = []

    #Get focal samples ( one for each simulation ) from input file 
    exp1 = str(gens) + 'gens' #Gets you the lines from bp file corresponding to the chosen pedigree based on gens
    exp2 = 'g' + str(gens+1) #Gets you the lines corresponding to the focal samples 

    bp_list_gens_all = flh.print_all_match(input_file, exp1, exp2) #List containing all strings corresponding to focal sample given gens
    bp_list_gens = bp_list_gens_all[sims[0]*2:sims[1]*2] #get sims corresponding to input
    #print(bp_list_gens[0],bp_list_gens[1])
    #Loop over the simulations
    for i,elem in enumerate(bp_list_gens):
        if i%2 != 0:
            continue

        print(i, len(bp_list_gens)) 
        #Ancestry parse the input strings corresponding to the current simulation. If simulation is weird and doesn't contain the right number of unique haplotypes, ancestry_parse returns an empty list and we continue
        focal_ancestry_bp = prs.ancestryParse([bp_list_gens[i],bp_list_gens[i+1]],admixture_type)
        if not focal_ancestry_bp:
            continue
        #Get dictionary of homozygous tracts for the current simulation 

        mapfile = '/fs/cbsubscb09/storage/resources/genetic_maps/refined_mf.simmap' #Reading in map_file to pass into haplo_parse
        map_df = pd.read_table(mapfile)

        focal_ancestry_parsed = [prs.haplo_parse(i,map_df) for i in focal_ancestry_bp]
        focal_ancestry_homoz_tracts = [prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'A'),\
                prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'B')] 

        #pull out tract lengths for each population, append to corresponding list, return lists 
        tract_lengths_A = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[0]] 
        tract_lengths_list_A.append(tract_lengths_A)
        tract_lengths_B = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[1]] 
        tract_lengths_list_B.append(tract_lengths_B)
    
    #if outfile is provided, pickle the lists for later use
    if outfile is not None:
        with open(outfile,"wb") as f:
            pickle.dump((tract_lengths_list_B,tract_lengths_list_B), f)

    #return lists 
    return(tract_lengths_list_A, tract_lengths_list_B)


def get_mle(input_list):
    '''
    Takes in list of lists of tract lenghts ( one for each simulation in bp file ) and returns mle estimate of admixture time for each sample as a list 
    '''
    mle_mean_list = []

    for sample_list in input_list:
        mle_mean_sample = sum(sample_list)/(2*len(sample_list))
        mle_mean_list.append(mle_mean_sample)    
    
    return(mle_mean_list)


if __name__ == "__main__":
    

