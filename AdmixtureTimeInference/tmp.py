
def extract_homozygous_tracts(hap1_dicts, hap2_dicts, haplotype_id):
    '''
    Takes in two haplo_parsed lists of dictionaries and the haplotype id (e.g 'CEU' or 'A' ) and returns list of dictionaries of ROHs for the chosen haplotype
    '''
    
#{{{
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

#}}}        
    #Return list of dicts containing homozygous tracts for input haplotype
    return(overlap_dicts)


def physical_to_genetic(bp_list, mapfile):
    '''
    Takes in input bp file as a list of strings and mapfile location and returns strings with genetic instead of physical positions 
    '''
#{{{
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
#}}}
    return(output_list)



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


