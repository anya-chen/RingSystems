import os
from save_ring_systems import save_ring_systems
from get_NP_rings import get_ring_systems_of_molecule_file
import multiprocessing as mp
from datetime import datetime

'''
Script to obtain the ring systems for each file of the organism subsets. Uses the function get_ring_systems_of_molecule_file of the get_NP_rings file. 
Use after assign_ids_organismSets.py
'''

def get_ring_systems(inputDir, file, OutputDir):
    '''
    Receives the ring systems for each file in the directory and saves them using save_ring_systems.

    :param dir: Directory with preprocessed molecule files
    :param file: a single tab separated file from the input directory with the columns
     index (should have not saved this), preprocessedSmiles, id, conID, mol_stereoIDs
    :param OutputDir: Directory for the result files. The absolute path to the single file is created according to the
    file name.
    '''
    ringsDF = get_ring_systems_of_molecule_file(inputDir, file)
    save_ring_systems(ringsDF, OutputDir + '/' + file.split('.')[0] +'Rings')
    print('Finished processing file {}'.format(file))

    
    
if __name__ == '__main__':
    start_time = datetime.now()
    
    inputDir = "/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/smilesHasRing_ids_optimized"
    outputDir = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/ringSystems'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    jobList = []
    for file in os.listdir(inputDir):
        jobList.append((inputDir, file, outputDir))

    pool = mp.Pool(processes=41)
    pool.starmap(get_ring_systems, jobList)
    pool.close()
    pool.join()
    
    print('Finished in:')
    print(datetime.now()-start_time)
