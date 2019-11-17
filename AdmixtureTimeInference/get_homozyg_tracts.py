"""This is the main module for AdmixtureTimeInference project

This module contains two functions: get_tract_lengths finds homozygosity tracts in simulated data and outputs tract lengths, while get_mle inference time since admixture from those tract lengths
"""

def get_tract_lengths(input_file, admixture_type, drop_ends=False):
    '''
    Args:
    input_file : Takes in an input .bpgen file (this is a .bp file with genetic instead of physical positions) 
    admixture type: list of founder populations (e.g ['0' , '1' , '0' , '1' ] ) 
    drop_ends : allows you to specify if segment at end of chromosome is to be included

    Returns:  
    returnDict: Dictionary with population identifiers as keys continaing the lists of lists. (One list for each simulation in input .bpgen file)
    returnDict: If mode = 'all' returnDict contains all tracts with keys 'hom_0','hom_1','het'
    '''

    print('Parameters printed from main.py: input_file = {}, admixture_type = {}, drop_ends = {}'.format(input_file, admixture_type, drop_ends))
            
    with open(input_file) as f:
        infile_lines = f.readlines()

    gens = int(re.search('(\d)gens', infile_lines[0]).group(1))
    assert 2**(gens) == len(admixture_type) 

    #tract_lengths_list_A = [] #initialize tract lengths list
    #tract_lengths_list_B = []
    tract_lengths_all_dicts = []

    #Loop over the simulations in inputfile
    for i in range(len(infile_lines)):
#{{{ 
        if i%2 != 0: #we want to loop in batches of two to process two haplotypes together
            continue

        #Ancestry parse the input strings corresponding to the current simulation.        
        
        try: 
            focal_ancestry_bp = prs.ancestryParse([infile_lines[i],infile_lines[i+1]],admixture_type)
        except:
            print('Error: Check that input file has even number of lines')
            exit()

       
       if drop_ends == True:
            focal_ancestry_parsed = [prs.haplo_parse(i,drop_ends=True) for i in focal_ancestry_bp]
        else:
            focal_ancestry_parsed = [prs.haplo_parse(i) for i in focal_ancestry_bp]
        
        
        #Get dictionary of homozygous tracts for the current simulation (if mode = 'all' get dictionary contiaining all tracts in different format)
        #UPDATE:mode option dropped. Compute all tracts by default !!!!!!

        print('extracting all tracts')
        
        focal_ancestry_all_tracts = prs.extract_all_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1])
        tract_lengths_all_dicts.append(focal_ancestry_all_tracts)
        
        #else:
        #    print('get_tract_lengths is running in mode "homzyg"; extract homozygous tracts only')
        #    focal_ancestry_homoz_tracts = [prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], '0'),\
        #    prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], '1')] 

        #    #pull out tract lengths for each population, append to corresponding list, return lists 
        #    tract_lengths_A = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[0]] 
        #    tract_lengths_list_A.append(tract_lengths_A)
        #    tract_lengths_B = [dicT['length'] for dicT in focal_ancestry_homoz_tracts[1]] 
        #    tract_lengths_list_B.append(tract_lengths_B)

#}}}

    #tract_lengths_homzyg_dicts = {"0":tract_lengths_list_A,"1":tract_lengths_list_B}
    returnDict = tract_lengths_all_dicts 

    return(returnDict)


if __name__ == "__main__":
#{{{ 



    import file_handling as flh
    import parsers as prs
    import numpy as np 
    import pandas as pd
    import argparse
    import re 
    import ast 

    # initiate the parser
    parser = argparse.ArgumentParser()

    # add long and short argument
    parser.add_argument("--inputfile", "-i",required=True, help="set input file (must be a bp file)")
    parser.add_argument("--admixture_label","-a",required=True, help='admixture type as key to admixtrue_labels dicFionary')
    parser.add_argument("--labels_file", "-l", required=True,help="admixture_labels as a file containing a dictionary string to evalutated literally")
    parser.add_argument("--drop_ends","-d", default=False, help='Specify as True if last blocks are to be dropped')
    parser.add_argument("--outfile", "-o", required=True,help="e.g outfile.pkl")
    parser.add_argument("--mode", "-m", default='all',required=True,help="specify whether to run in mode 'all' or 'homzyg' ")
    # read arguments from the command line
    args = parser.parse_args()



    #~~~~~~~~ Body ~~~~~~~~~~~~~~~~~ #
    #Read in labels file literally (contains dictionary)
    with open(args.labels_file) as f:
        dict_string = f.readline()

    label_dict = ast.literal_eval(dict_string)

    admixture_type_string = label_dict[str(args.admixture_label)] #Read in admixture type string from input dictionary give input admixture type label
    admixture_type_list = list(admixture_type_string)
    
    #Strip header column from infile     
    header_list=[re.findall(r'(.*?\s)1\|',line)[0] for line in open(args.inputfile)]   

    #Final tract length options 
    if args.mode=='homzyg':
        res = get_tract_lengths(args.inputfile, admixture_type_list, args.drop_ends, mode=args.mode)

        #Write to outfile (one column for A, one for B) with header column 
        with open(args.outfile,mode='w') as f:
            f.write("header\thom_0\thom_1\n")
            for indx, elem in enumerate(res['0']):
                f.write("{}\t{}\t{}\n".format(header_list[indx*2],res['0'][indx],res['1'][indx])) 

    if args.mode=='all':
        res = get_tract_lengths(args.inputfile, admixture_type_list, args.drop_ends, mode=args.mode)

        #Write to outfile (one column for A, one for B)  with header column 
        with open(args.outfile,mode='w') as f:
            f.write("header\thom_0\thom_1\thet\thom0_status\thom1_status\thet_status\n")
            for indx, dicT in enumerate(res):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header_list[indx*2],dicT['hom_0'], dicT['hom_1'], dicT['het'], dicT['hom0_status'], dicT['hom1_status'], dicT['het_status'])) 


#}}}
