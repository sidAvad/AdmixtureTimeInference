import parsers as prs 
import importlib

#Testing parsers module 
importlib.reload(prs)
input_strings = prs.print_two_matches('/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Output/ped-sim/admixture-time.bp', "3gens","g4")

input_strings_ancestry = prs.ancestryParse(input_strings, ['A','B','A','B','A','B','A','B'])

input_dicts = [None]*2
input_dicts[0] = prs.haplo_parse(input_strings_ancestry[0])
input_dicts[1] = prs.haplo_parse(input_strings_ancestry[1])

importlib.reload(prs)
prs.extract_homozygous_tracts(input_dicts[0], input_dicts[1], 'B')

#Testing main module 
import main

importlib.reload(main)
input_file = '/fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/Output/ped-sim/admixture-time.bp'
res = main.get_tract_lengths(input_file, 4, ['A','B','A','B','A','B','A','B','A','B','A','B','A','B','A','B'])



#Testing haploparse function with genetic_map option
importlib.reload(main)
importlib.reload(prs)
mapfile ='/fs/cbsubscb09/storage/resources/genetic_maps/refined_mf.simmap'
map_df = pd.read_table(mapfile)
input_strings_ancestry = prs.ancestryParse(input_strings, ['A','B','A','B','A','B','A','B'])

output_dicts = [None]*2
output_dicts[0] = prs.haplo_parse(input_strings_ancestry[0], map_df)
output_dicts[1] = prs.haplo_parse(input_strings_ancestry[1], map_df)

importlib.reload(prs)
res = main.get_tract_lengths(input_file,4,['A','B']*8,[0,10],'results_4gens_batch1.pkl')
mle_t = main.get_mle(res[0])
