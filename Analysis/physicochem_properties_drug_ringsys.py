from physicochem_properties import *

inputfile = '/data/local/ringsys/202111_approveddrug/uniqueRingSystems/drug_uniqueRingSystems_noStereo.smi'
outputdir = '/data/local/ringsys/20211202_analysis/'

ringsDF = get_dataframe_from_ring_csv(inputfile)

ringsDF.loc[:, 'db'] = 'Drugs'

ringsDF = get_physicochemical_properties(ringsDF)

largest_ringsystem_by_column(ringsDF, 'heavy_atoms', outputdir)
largest_ringsystem_by_column(ringsDF, 'numRings', outputdir)

save_rings_df_with_properties(ringsDF, outputdir + 'drugs_all_descriptors.csv')
calculate_mean_and_sd(ringsDF, outputdir + 'drugs_mean_sd_results.txt')
print('Finished calculating properties for drugs ring systems')
