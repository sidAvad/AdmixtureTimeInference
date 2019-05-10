#Importing libraries
import file_handling as flh
import parsers as prs 
import importlib
import main


#arguments
meios = 2
input_file = '/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Output/ped-sim/admixture-time_' + str(meios) + 'meios_poisson_sexavg.bp'

exp1 = 'g' + str(meios+2) #Gets you the lines corresponding to the focal 
exp2 =  'meios288' 
bp_list_focal = flh.print_all_match(input_file, exp1 ,exp2)
bp_list_focal = flh.print_two_matches(input_file, exp1,exp2) #List containing all strings corresponding to focal sample given gens


mapfile ='/fs/cbsubscb09/storage/resources/genetic_maps/refined_mf.simmap' 
importlib.reload(prs)

#res = prs.phys2gen_string(bp_list_focal[0],mapfile)
res_list = prs.phys2gen_list(bp_list_focal,mapfile)

admixture_type = 'ABABABAB'

importlib.reload(prs)
focal_ancestry_bp = prs.ancestryParse(res_list,admixture_type)
focal_ancestry_parsed = [prs.haplo_parse(i) for i in focal_ancestry_bp]
focal_ancestry_parsed_cut = [prs.haplo_parse(i,drop_ends=True) for i in focal_ancestry_bp]

focal_ancestry_homoz_tracts = [prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'A'),\
                               prs.extract_homozygous_tracts(focal_ancestry_parsed[0], focal_ancestry_parsed[1], 'B')] 

importlib.reload(prs)
collapsed = prs.collapse_blocks(focal_ancestry_parsed[0])

#######################
input_strings_ancestry = prs.ancestryParse(input_strings, ['A','B','A','B','A','B','A','B'])

meios = 1
input_file = '/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Output/ped-sim/admixture-time_' + str(meios) + 'meios_poisson_sexavg.bp'
mapfile ='/fs/cbsubscb09/storage/resources/genetic_maps/refined_mf.simmap' 

importlib.reload(prs)
importlib.reload(main)
res = main.get_tract_lengths(input_file,mapfile,['A','B']*2 , [0,2])
mle_t = main.get_mle(res[0])

######################

import pandas as pd 
import re
import parsers2 as prs 

mapfile ='/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Input/averaged_mf.simmap' 
bpfile = '/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Output/ped-sim/admixture-time_1meios_intf_sexspf.bp'
bpfile_test = '/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/bpfile_test.bp'

importlib.reload(prs)
res = prs.phys2gen_file(bpfile_test, mapfile)


