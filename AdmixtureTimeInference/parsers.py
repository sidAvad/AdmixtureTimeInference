#functions found in Tools/bpTools
#Import required modules
import re   
import pandas as pd 

def find_genetic_position(dataframe, physical_position):
    '''
     Takes in genetic_map dataframe ( containing sex_averaged column ) and physical_position and returns corresponding interpolated genetic position
    '''
    df = dataframe.copy(deep=True) 

    if physical_position <= df['pos'].iloc[0]:
        genetic_position = round(float(df['sex_averaged'].iloc[0]),6)
        return(genetic_position)

    if physical_position >= df['pos'].iloc[-1]:
        genetic_position = round(float(df['sex_averaged'].iloc[-1]),6)
        return(genetic_position)

    df.loc[-1, 'pos'] = physical_position
    df  = df.sort_values('pos').reset_index(drop=True)
    df = df.interpolate()
    genetic_position = round(float(df['sex_averaged'][df.pos == physical_position].iloc[0]),6) #sometimes physical position matches more than one location and the resulting series cannot be converted to a float 
      
    return(genetic_position)



def haplo_parse(input_string , map_df=False):
    '''
    Args:
    input_string : Takes input string in bp format (output from ped-sim) 
    genetic_map : Replaces physical positions with genetic position based on input map dataframe

    Returns:
    haplo_dicts : list of dictionaries containing chromosome number, haplotype,
                  block length and [start, stop] positions as keys 

    '''
    
    input_string = input_string.replace('\n',' ') #Replace trailing newline with a space for easier regex processing
    parsed_list = re.findall(r'(\d{1,2}\|\d+\s(?:[A-Za-z0-9]+:\d+\s)+)', input_string, re.MULTILINE) #Killer regex here. 
    haplo_lists_of_dicts = [None]*len(parsed_list)

    if not isinstance(map_df, pd.DataFrame): #if map_df not provided, we use physical positions ( default value of map_df is false )

        for i in range(len(parsed_list)): #Loop over the chromosomes
            #regex parse each list entry to pull out all the information we need 
            chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
            chrom_start = int(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
            haplotypes = re.findall(r'\s(.*?):', parsed_list[i]) 
            chrom_end_list = re.findall(r':(.*?)\s', parsed_list[i])
            range_list = [[int(chrom_start),int(chrom_end_list[i])]  if i==0 else [int(chrom_end_list[i-1]),int(chrom_end_list[i])] for i in range(len(chrom_end_list))]
            chrom_list = [{'haplotypes':haplotypes[j] ,'startstop':range_list[j], 'chromosome':chrom, 'length':range_list[j][1] - range_list[j][0] } for j in range(len(haplotypes)) ] 
            haplo_lists_of_dicts[i] = chrom_list # Assign chromosome i's list of haplotype dictionaries to a list 

    else: #If map_df is provided, we replace phsyical positions with genetic

        for i in range(len(parsed_list)): #Loop over the chromosomes
            #regex parse each list entry to pull out all the information we need 
            chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
            chrom_start = int(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
            
            #Get df correpsonding to current chromosome 
            chrom_df = map_df[map_df['#chr'] == int(chrom)]
            chrom_df.loc[:,'sex_averaged'] = (chrom_df.male_cM + chrom_df.female_cM)/2

            #convert chrom_start into genetic position
            chrom_start_genetic = find_genetic_position(chrom_df, chrom_start)

            #convert chrom_end_list into genetic positions
            chrom_end_list = re.findall(r':(.*?)\s', parsed_list[i])
            chrom_end_genetic_list = [find_genetic_position(chrom_df, int(i)) for i  in chrom_end_list]
            range_list = [[chrom_start_genetic,chrom_end_genetic_list[i]]  if i==0 else [chrom_end_genetic_list[i-1],chrom_end_genetic_list[i]] for i in range(len(chrom_end_genetic_list))]

            haplotypes = re.findall(r'\s(.*?):', parsed_list[i]) 
            chrom_list = [{'haplotypes':haplotypes[j] ,'startstop':range_list[j], 'chromosome':chrom, 'length':range_list[j][1] - range_list[j][0] } for j in range(len(haplotypes)) ] 
            #print(chrom_list)
            haplo_lists_of_dicts[i] = chrom_list # Assign chromosome i's list of haplotype dictionaries to a list 


    haplo_dicts = [item for sublist in haplo_lists_of_dicts for item in sublist] #Collapse list 
    return(haplo_dicts) 



def ancestryParse(input_strings, admixed_branches):
    '''
    Args:
    input_strings : Takes input strings - one for each chromosome - in bp format (output from ped-sim) 
    admixed_branches : list of admixture types for each branch ( list length should be equal to number of founders)

    Returns:
    output_strings : list of strings (one element for each chromosome) in bp format with haplotype id replaced by population/ancestry ids 
    '''

    #get the focal samples from `gens`  
    #focal_gen = 'g' + str(gens + 1)  

    #Get all haplotypes in the input bp strings 
    haplotypes = [re.findall(r'\s(\d+):' , i) for i in input_strings]
    haplotypes_unique_ordered = list(set([int(y) for x in haplotypes for y in x ]))
    #Check if admixed_branches is compatible with gens else throw exception and return an empty list
    try: 
        assert 2*len(admixed_branches) == len(haplotypes_unique_ordered)
    except:
        print("Warning: Incorrect number of unique haplotypes found given admixed_branches, returning empty list")
        return([])

    #Replace input string haplotypes with population identifiers based on admixed_branches
    #Create dictionary mapping haplotypes to admixture types 
    admixed_branches_expanded = [admixed_branches[i//2] for i in range(len(admixed_branches)*2)] #Duplicate input for easier mapping to haplotypes in bp file 
    
    #Replace input strings with corresponding admixture type
    output_strings = [None]*2
    for j in range(len(input_strings)): 
        output_strings[j] = input_strings[j]
        for i in range(len(haplotypes_unique_ordered)): 
            output_strings[j] = output_strings[j].replace(str(haplotypes_unique_ordered[i]) + ':', admixed_branches_expanded[i] + ':')

    #input_string = input_string.replace('\n',' ') #Replace trailing newline with a space for easier regex processing
    return(output_strings)


def extract_homozygous_tracts(hap1_dicts, hap2_dicts, haplotype_id):
    '''
    Takes in two haplo_parsed lists of dictionaries and the haplotype id (e.g 'CEU' or 'A' ) and returns list of dictionaries of ROHs for the chosen haplotype
    ''' 

    #Initialize overlap_dicts
    output_dicts = [] #Initialize return list 
    #Loop over chromosomes 
    for chromosome in range(1,23): 

        #Subset for chromosome and haplotype_id
        hap1_chr = [ dicT for dicT in hap1_dicts if dicT['chromosome'] == str(chromosome) and dicT['haplotypes'] == haplotype_id ]
        hap2_chr = [ dicT for dicT in hap2_dicts if dicT['chromosome'] == str(chromosome) and dicT['haplotypes'] == haplotype_id ]
        #Loop over haplotype blocks for hap1 in the current chromosome
        for i, dictionary1 in enumerate(hap1_chr):
            #if no more hap1 blocks left break and move to next chromosome
            if i == len(hap1_chr):
                break
            #Loop over haplotype_ids in hap2 for current chromosome
            for j, dictionary2 in enumerate(hap2_chr):
            #if no more hap2 blocks left, break and move outer loop 
                if j == len(hap2_chr):
                    break
                #Get start and stop positions of blocks in each haplotype
                hap1_start = dictionary1['startstop'][0]
                hap1_end = dictionary1['startstop'][1]
                hap2_start = dictionary2['startstop'][0] 
                hap2_end = dictionary2['startstop'][1]
                #if hap2 is completely behind hap1, move inner loop
                if hap2_end < hap1_start:
                    #print('2 behind 1:',dictionary1['startstop'], dictionary2['startstop'])
                    continue
                #else if hap2 is completely ahead of hap1, break and move outer loop 
                elif hap1_end < hap2_start:
                    #print('1 behind 2:',dictionary1['startstop'], dictionary2['startstop'])
                    break
                #if overhang in hap2, compute overlap dict and move outer loop
                elif hap2_end > hap1_end:
                    #Compute overlap_dict and append to output_list 
                    overlap_dict = {'haplotypes':None, 'overlap':None, 'chromosome':None, 'length':None} #initialize overlap_dict
                    overlap_dict['haplotypes'] = dictionary1['haplotypes'] + dictionary2['haplotypes']  
                    overlap_dict['overlap'] = [max(hap1_start,hap2_start), min(hap1_end,hap2_end)]
                    overlap_dict['chromosome'] = dictionary1['chromosome']
                    overlap_dict['length'] = min(hap1_end,hap2_end) - max(hap1_start,hap2_start)
                    output_dicts.append(overlap_dict)
                    #print(dictionary1['startstop'], dictionary2['startstop']) #Print statements for debugging
                    #print(overlap_dict)
                    break
                #else , compute overlap dict and move inner loop 
                else: 
                    overlap_dict = {'haplotypes':None, 'overlap':None, 'chromosome':None, 'length':None} #initialize overlap_dict
                    overlap_dict['haplotypes'] = dictionary1['haplotypes'] + dictionary2['haplotypes']  
                    overlap_dict['overlap'] = [max(hap1_start,hap2_start), min(hap1_end,hap2_end)]
                    overlap_dict['chromosome'] = dictionary1['chromosome']
                    overlap_dict['length'] = min(hap1_end,hap2_end) - max(hap1_start,hap2_start)
                    output_dicts.append(overlap_dict)
                    #print(dictionary1['startstop'], dictionary2['startstop'])
                    #print(overlap_dict)
                    continue

    #Return list of dicts containing homozygous tracts for input haplotype
    return(output_dicts)




#DATA WRANGLING FUNCTIONS
def print_two_matches(input_file, expression1, expression2):
    '''
    return  list of strings in input bp file containing the first two matches to both input expressions
    '''
    infile = open(input_file,'r') 
    num_found = 0 
    found_matches = [None]*2

    for line in infile:
        if re.search(expression1, line):
            if re.search(expression2, line): #Not the most efficient way to 'and' two regexs
                found_matches[num_found] = line.strip('\n')
                num_found = num_found + 1  #break after matching both chromosomes
                if num_found == 2: 
                    break

    infile.close()
    return(found_matches)

def print_all_match(input_file, expression1, expression2):
    '''
    return list with all lines containing matches to expression1 and expression2 in input_file
    '''
    infile = open(input_file,'r') 

    match_list = list()

    for line in infile:
        if re.search(expression1 , line):
            if re.search(expression2, line): #Not the most efficient way to 'and' two regexs
                match_list.append(line.strip('\n'))

    infile.close()
    return(match_list)


