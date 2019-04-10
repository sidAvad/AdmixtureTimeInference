
def extract_homozygous_tracts(hap1_dicts, hap2_dicts, haplotype_id):
    '''
    Takes in two haplo_parsed lists of dictionaries and the haplotype id (e.g 'CEU' or 'A' ) and returns list of dictionaries of ROHs for the chosen haplotype
    '''
    
    start = 0
    overlap_dicts = [] #Initialize return list 
    #Loop over all the chosen haplotypes in chrom1
    for dictionary in hap1_dicts:

        if dictionary['haplotypes'] is not haplotype_id: #Check that use specified input is the one being tested
            continue

        #While loop over chosen haplotypes in hap2 
        position = start #Pick up from last unmatched haplotype in hap2 
        stop = False  #Get the while loop going
        
        while not stop:
            if dictionary['haplotypes'] is not hap2_dicts[position]['haplotypes']: #Check that the same haplotype id is being compared
               position = position + 1
               continue
            
            print(dictionary['startstop'], dictionary['haplotypes'])
            print(position)
            #Find overlap, save to ROH dict and move to next haplotype
            #Find overlap 
            r1_start = dictionary['startstop'][0]
            r1_end = dictionary['startstop'][1]
            r2_start = hap2_dicts[position]['startstop'][0] 
            r2_end = hap2_dicts[position]['startstop'][1] 

            if r2_end <  r1_start: #Continue on hap2if no overlap and hap2 range behind hap1 range 
                position = position + 1
                continue

            if r2_start > r1_end: #When segment in hap2 is beyond hap1 segment, it is time to move the outer loop
                start = position
                break

            overlap_dict = {'haplotypes':None, 'overlap':None, 'chromosome':None, 'length':None} #initialize overlap_dict

            #update output list of dicts with homozygous tract information dict 
            overlap_dict['haplotypes'] = dictionary['haplotypes'] + hap2_dicts[position]['haplotypes']  
            overlap_dict['overlap'] = [max(r1_start,r2_start), min(r1_end,r2_end)]
            overlap_dict['chromosome'] = dictionary['chromosome']
            overlap_dict['length'] = min(r1_end,r2_end) - max(r1_start,r2_start)

            print(overlap_dict) 
            overlap_dicts.append(overlap_dict)

            #End while loop if haplotype is partially overlapping at end or completely non-overlapping at the end , set stop to TRUE and set start to current position for the next loop 
            if r2_end >= r1_end:
                stop = True #End the while loop 
                start = position
            else:
                position = position + 1

        
    #Return list of dicts containing homozygous tracts for input haplotype
    return(overlap_dicts)


def physical_to_genetic(bp_list, mapfile):
    '''
    Takes in input bp file as a list of strings and mapfile location and returns strings with genetic instead of physical positions 
    '''
    #Read in map file 
    genetic_map = pd.read_table(mapfile)    
    genetic_map['sex_averaged'] = (genetic_map.male_cM + genetic_map.female_cM)/2 #Get sex averaged genetic positions
    
    #Conver input_list physical positions to genetic positions and store in new list
    output_list = []
    for string in bp_list:
        string = string.replace('\n',' ') #Replace trailing newline with a space for easier regex processing
        parsed_list = re.findall(r'(\d{1,2}\|\d+\s(?:[A-Za-z0-9]+:\d+\s)+)', string, re.MULTILINE) #Killer regex here. 
        output_list.append(parsed_list)

    #return list of genetic positions

    return(output_list)
