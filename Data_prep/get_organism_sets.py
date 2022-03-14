import csv
from pathlib import Path
'''
This script is used to get sets of NPs from different organisms from the refined coconut.

'''

refinedCoconut = '/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/refinedCoconut.sourceNP.csv'
outputDir = Path(refinedCoconut).parent

plants_sources = ['tppt', 'cmaup', 'respect', 'hitdb', 'biofacquim', 'etmdb', 'mitishamba','vietherb', 'afromalariadb', 'afrodb', 'npact', 'himdb', 'afrocancer', 'p-anapl', 'spektraris', 'knapsack', 'inpacdb', 'tipdb','tcmid']
bacteria_sources = ['np_atlas_2019_12', 'piellabdata', 'streptomedb3', 'knapsack', 'npatlas', 'streptomedb']
fungi_sources = ['np_atlas_2019_12', 'biofacquim', 'npatlas', 'knapsack']
marine_sources = ['mnp', 'cmnpd']

plants_set = outputDir / 'plants.csv'
bacteria_set = outputDir / 'bacteria.csv'
fungi_set = outputDir / 'fungi.csv'
marine_set = outputDir / 'marine.csv'

with open(plants_set, 'w') as plants:
    with open(bacteria_set, 'w') as bacteria:
        with open(fungi_set, 'w') as fungi:
            with open(marine_set, 'w') as marine:
                with open(refinedCoconut, 'r') as f:
                    moleculeCsv = csv.reader(f, delimiter='\t')
                    header = next(moleculeCsv)
                    print(header)
                    plants.write( '\t'.join(map(str, header)))
                    plants.write('\n')
                    bacteria.write( '\t'.join(map(str, header)))
                    bacteria.write('\n')
                    fungi.write( '\t'.join(map(str, header)))
                    fungi.write('\n')                    
                    marine.write( '\t'.join(map(str, header)))
                    marine.write('\n')

                    for row in moleculeCsv:
                        if row[1] in plants_sources and ('"plants"' in row[5]):
                            plants.write('\t'.join(map(str,row)))
                            plants.write('\n')
                        if row[1] in bacteria_sources and ('"bacteria"' in row[5]):
                            bacteria.write('\t'.join(map(str,row)))
                            bacteria.write('\n')
                        if row[1] in fungi_sources and ('"fungi"' in row[5]):
                            fungi.write('\t'.join(map(str,row)))
                            fungi.write('\n')
                        if row[1] in marine_sources:
                            marine.write('\t'.join(map(str,row)))
                            marine.write('\n')
            marine.close()
        fungi.close()
    bacteria.close()
plants.close()
