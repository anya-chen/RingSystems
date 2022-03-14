import os
from save_ring_systems import save_ring_systems
from get_NP_rings import get_ring_systems_of_molecule_file
import multiprocessing as mp
from datetime import datetime

'''
Script to obtain the ring systems for each file of the ZINC database. Uses the function get_ring_systems_of_molecule_file
of the get_NP_rings file. As the number of all ring systems would exceed the memory, the unique ring system set is
obtained using get_unique_ZINC_ring_system_set.py.
Use after assign_ids_molecules.py
'''


def get_ring_systems(dir, file, OutputDir):
    '''
    Receives the ring systems for each file in the directory and saves them using save_ring_systems.

    :param dir: Directory with preprocessed molecule files
    :param file: a single tab separated file from the input directory with the columns
     index (should have not saved this), preprocessedSmiles, id, conID, mol_stereoIDs
    :param OutputDir: Directory for the result files. The absolute path to the single file is created according to the
    file name.
    '''
    ringsDF = get_ring_systems_of_molecule_file(dir, file)
    save_ring_systems(ringsDF, OutputDir + '/Rings_' + file.split('.')[0])
    print('Finished processing file {}'.format(file))


if __name__ == '__main__':
    start_time = datetime.now()
    
    ZINCdir = "/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/smilesHasRing_ids_optimized"
    ZINCoutput = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/ringSystems'
    if not os.path.exists(ZINCoutput):
        os.makedirs(ZINCoutput)

    jobList = []
    for file in os.listdir(ZINCdir):
        if "Rings_" + file not in os.listdir(ZINCoutput):
            jobList.append((ZINCdir, file, ZINCoutput))

    pool = mp.Pool(processes=14)
    pool.starmap(get_ring_systems, jobList)
    pool.close()
    pool.join()
    
    print('Finished in:')
    print(datetime.now()-start_time)
