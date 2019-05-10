#Import required modules
import re   
import pandas as pd 
import numpy as np 


def haplo_parse(input_string, collapse = True, drop_ends=False):
    '''
    Args:
    input_string : Takes input string in bp format (output from ped-sim) 
    collapse : Collapses identical adjacent blocks by default
    drop_ends : drops final segment to avoid chromosome end bias. False by default
    Returns:
    haplo_dicts : list of dictionaries containing chromosome number, haplotype,
                  block length and [start, stop] positions as keys 

    '''
 #{{{ 
    input_string = input_string + ' ' #Add space for better regex processing 
    parsed_list = re.findall(r'\d+\|.+?(?=\d{1,2}\||$)', input_string, re.MULTILINE) #Killer regex here. Captures information chromosome by chromosome 
    haplo_lists_of_dicts = [None]*len(parsed_list)

    for i in range(len(parsed_list)): #Loop over the chromosomes
        #regex parse each list entry to pull out all the information we need 
        chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
        chrom_start = float(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
        haplotypes = re.findall(r'\s(.*?):', parsed_list[i]) 
        chrom_end_list = re.findall(r':(.*?)\s', parsed_list[i])
        range_list = [[float(chrom_start),float(chrom_end_list[i])]  if i==0 else [float(chrom_end_list[i-1]),float(chrom_end_list[i])] for i in range(len(chrom_end_list))]

        chrom_list = [{'haplotypes':haplotypes[j] ,'startstop':range_list[j], 'chromosome':chrom, 'length':range_list[j][1] - range_list[j][0] } for j in range(len(haplotypes)) ] 

        haplo_lists_of_dicts[i] = chrom_list # Assing chromosome i's list of haplotype dictionaries to a list 

    haplo_dicts = [item for sublist in haplo_lists_of_dicts for item in sublist] #Collapse list 
   
    #If option collapse is provided ( provided by default) we collapse contiguous blocks and return list 
    if collapse:
       haplo_dicts = collapse_blocks(haplo_dicts)

    #If option drop_ends is provided ( not provided by default), drop final block for each chromosome
    if drop_ends==True: 

        print('printed from haplo_parse - Dropping ends')
        haplo_dicts_drop_lists = []
        for chr in range(1,23): #loop over chromosomes
            haplo_dicts_chr = [elem for elem in haplo_dicts if elem['chromosome'] == str(chr)]
            del haplo_dicts_chr[-1]
            haplo_dicts_drop_lists.append(haplo_dicts_chr) 

        haplo_dicts_drop = [item for sublist in haplo_dicts_drop_lists for item in sublist]
        return(haplo_dicts_drop)

        # }}}
    return(haplo_dicts) 


def collapse_blocks(haplodicts):
    '''
    Collapses contiguous block in haplo-parsed dictionary
    '''
#{{{
    #Loop over blocks in haplodicts
    contig_list = []
    extending = False
    for index,elem in enumerate(haplodicts):  
    
        #While extending, if next haplotype is not the same as current OR not on the same chromosome, end extended block at current, else go to next block
        if extending:
            #First we check if we've reached the end 
            if index == len(haplodicts) - 1: #If we've extended into last block, save end state and continue
                contig_block_end = elem['startstop'][1]

                new_elem = {'haplotypes':contig_block_haplotype,'startstop':[contig_block_start,contig_block_end], 'chromosome':contig_block_chromosome, 'length':contig_block_end - contig_block_start}
                contig_list.append(new_elem)

                extending = False
                continue

            if elem['haplotypes'] != haplodicts[index+1]['haplotypes'] or elem['chromosome'] != haplodicts[index+1]['chromosome'] : #If we've extended into different block (different chromosome or different haplotype ), save end state and continue 
                contig_block_end = elem['startstop'][1]

                new_elem = {'haplotypes':contig_block_haplotype,'startstop':[contig_block_start,contig_block_end], 'chromosome':contig_block_chromosome, 'length':contig_block_end - contig_block_start}
                contig_list.append(new_elem)

                extending = False
                continue

            else:
                continue

        #If we've reached the end, and are not extending, append current element to final list 
        if index == len(haplodicts) - 1:
            contig_list.append(elem)
            continue

        #If next block is same as current, AND on the same chromosome :  Save start, chromosome and haplotype information and start extending. Otherwise append current element to final list and continue  
        if elem['haplotypes'] == haplodicts[index+1]['haplotypes'] and elem['chromosome'] == haplodicts[index+1]['chromosome']:     
            contig_block_start = elem['startstop'][0]
            contig_block_haplotype = elem['haplotypes'] 
            contig_block_chromosome = elem['chromosome'] 
            extending = True
        else:
            contig_list.append(elem)
            continue
            #}}}
    return(contig_list)

def ancestryParse(input_strings, admixed_branches):
    '''
    Args:
    input_strings : Takes input strings - one for each chromosome - in bp format (output from ped-sim) 
    admixed_branches : list of admixture types for each branch ( list length should be equal to number of founders)

    Returns:
    output_strings : list of strings (one element for each chromosome) in bp format with haplotype id replaced by population/ancestry ids 
    '''
#{{{

    #Get founder haplotypes using modular arithmetic ( based on number of meises simulated and current simulation number ) 
    sim_no = int(re.search('(\d+)_',input_strings[0]).group(1))
    num_founders = 2**(int(input_strings[0][0]) + 1)
    num_haplotypes = 2*num_founders

    assert len(admixed_branches) == num_founders # Check that both input strings correspond to admixed_branches input by checking how many meoises simulated ( 1st character of bp file will give us this information ) 


    first_id = (sim_no - 1)*num_haplotypes
    last_id = sim_no*num_haplotypes - 1 
    
    founder_haplotypes = [str(x) for x in range(first_id, last_id+1)]

    #Replace input string haplotypes with population identifiers based on admixed_branches
    ##Create dictionary mapping haplotypes to admixture types 
    admixed_branches_expanded = [admixed_branches[i//2] for i in range(len(admixed_branches)*2)] #Duplicate input for mapping to haplotypes in bp file 
    
    ##Replace input strings with corresponding admixture type
    output_strings = [None]*2
    for j in range(len(input_strings)): 
        output_strings[j] = input_strings[j]
        for i in range(len(founder_haplotypes)): 
            output_strings[j] = output_strings[j].replace(' ' + str(founder_haplotypes[i]) + ':', ' ' + admixed_branches_expanded[i] + ':') #Extra space is to ensure exact matches ( e.g, so that 16: doesn't match 6: )
            #}}}
    return(output_strings)


def extract_homozygous_tracts(hap1_dicts, hap2_dicts, haplotype_id):
    '''Takes in two haplo_parsed lists of dictionaries and the haplotype id (e.g 'CEU' or 'A' ) and returns list of dictionaries of ROHs for the chosen haplotype
    '''
#{{{
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
#}}}
    #Return list of dicts containing homozygous tracts for input haplotype
    return(output_dicts)


def phys2gen_string(bpstring,mapfile_list):
    '''
    Converts bp string with physical positions to bp string with interpolated genetic positions given a mapfile read in as a list of strings
    '''

#{{{
    #Check if mapfile has appropriate number of columns
    try:
        len(mapfile_list[0]) == 4
    except Warning:
        print('mapfile does not have the expected number of columns ; please provide mapfile with sex-specific columns only') 

    #preliminary regex parsing of input string  
    bpstring = bpstring + ' '  #Replace trailing newline with a space for easier regex processing
    parsed_list = re.findall(r'(\d{1,2}\|\d+\s(?:[A-Za-z0-9]+:\d+\s)+)', bpstring, re.MULTILINE) #Killer regex here. 
    #output_list = []

    #Initialize lists with all final physical and genetic positions (from all chromosomes)
    physical_positions_all = []
    genetic_positions_all = [] 

    #Loop over chromosomes and get physical positions 
    for i in range(len(parsed_list)):
        chrom = re.search(r'^[^|]*',parsed_list[i]).group(0)
        chrom_start = int(re.search(r'\|(.*?)\s', parsed_list[i]).group(1))
        chrom_end_list = re.findall(r':(\d+)', parsed_list[i])
        #make chrom_start and chrom_end_list into one list and convert to ints 
        physical_positions = [chrom_start] + [int(x) for x in chrom_end_list]

        #Get mapfile_list entries corresponding to current chromosome
        mapfile_chr = [x for x in mapfile_list if x[0] == chrom]   
        #print(chrom, mapfile_chr[-1])

        #interpolate and get genetic position for each position
        xp = [int(x[1]) for x in mapfile_chr] 
        fp = [(float(x[2]) + float(x[3]))/2 for x in mapfile_chr] #Getting sex averaged genetic position entries
        genetic_positions = [np.interp(x, xp ,fp) for x in physical_positions] #Interpolating,rounding, and pulling out genetic positions for each phsyical position 
        #print(physical_positions, genetic_positions)
        #Append physical and genetic_positions for current chromosome to final list of physical and corresponding genetic positions 
        physical_positions_all = physical_positions_all + physical_positions
        genetic_positions_all = genetic_positions_all + genetic_positions


    #Checking that all elements of physical_positions_all are unique 
    try: 
        len(set(physical_positions_all)) == len(physical_positions_all) 
    except:
        print('repeating breakpoints across chromosomes. Returning list of physcial positions across all chromosomes. Check input')
        return(physical_positions_all)

    #Find and replace physical position with interpolated genetic position
    for match_obj in re.finditer(':\d+|\|\d+', bpstring):

        match = match_obj.group(0) #Get string from match object 

        if ':' in match:
            match_int = int(match.strip(':'))
            stripped_char = ':'
        if '|' in match:
            match_int = int(match.strip('|'))
            stripped_char = '|' 

        match_genpos = round(genetic_positions_all[physical_positions_all.index(match_int)],3) 
        replace_str = stripped_char + str(match_genpos)

        bpstring = bpstring.replace(match,replace_str, 1) 

#}}}
    #Return final bp string with trailing whitesapce removed
    return(bpstring[:-1])


def phys2gen_list(bplist,mapfile):
    '''Take in list of bpfile entries with physical positions and mapfile and returns list of bpfile entires with interpolated genetic positions
    '''

#{{{
    #read in mapfile as list of strings  
    with open(mapfile) as f:
       mapfile_list = f.readlines()

    mapfile_list = [x.strip().split('\t') for x in mapfile_list]

    outlist = [] 
    #Loop over input and call phys2gen_string to do the heavy lifting 
    for i,bpstring in enumerate(bplist):
        outstring = phys2gen_string(bpstring,mapfile_list)
        outlist.append(outstring)

#}}}
    return(outlist)

    
