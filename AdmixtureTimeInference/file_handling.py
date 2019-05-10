#Contains functions to find lines in a text file according to regex matches
import re 

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

#def print_all_match(input_file, expression1, expression2):
#    '''
#    return list with all lines containing matches to expression1 and expression2 in input_file
#    '''
#    infile = open(input_file,'r') 
#
#    match_list = list()
#
#    for line in infile:
#        if re.search(expression1 , line):
#            if re.search(expression2, line): #Not the most efficient way to 'and' two regexs
#                match_list.append(line.strip('\n'))
#
#    infile.close()
#    return(match_list)

def print_all_match(input_file, expression1):
    '''
    return list with all lines containing matches to expression1 in input_file
    '''

    infile = open(input_file,'r') 

    match_list = list()

    for line in infile:
        if re.search(expression1 , line):
           match_list.append(line.strip('\n'))

    infile.close()
    return(match_list)
