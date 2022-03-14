import os
from pathlib import Path
from rdkit import Chem
from get_ring_systems import process_molecules
from save_ring_systems import save_ring_systems
import pandas as pd
from datetime import datetime

'''
Script to obtain the ring systems of the natural products.
Use after assign_ids_molecules.py
'''

def get_NP_rings(NPdir, OutputDir):
    '''
    Takes all files with preprocessed molecules and saves their ring systems

    :param NPdir: directory with preprocessed natural product SMILES files
    :param OutputDir: directory to save the ring systems of each file
    '''
    for file in os.listdir(NPdir):
        print(file)
        filename = file.split('.')[0][:-11]
        ringsDF = get_ring_systems_of_molecule_file(NPdir, file)
        save_ring_systems(ringsDF, str(OutputDir) + "/" + filename + "Rings")


def get_ring_systems_of_molecule_file(dir, file):
    '''
    Function to obtain the ring systems of the molecules of a single file

    :param dir:Directory of molecule file
    :param file: file with preprocessed molecules obtained by preprocess_molecules.py
    :return:dataframe with ring systems of inputmolecules
    '''
    filepath = os.path.join(dir, file)
    
    ringsDF = process_molecules(filepath)
    return ringsDF



if __name__ == "__main__":   
    start_time = datetime.now()

    inputPath = "/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/smilesHasRing_ids_optimized"
    saveRingsPath = Path(inputPath).parent / 'ringSystems'
    if not saveRingsPath.exists():
        saveRingsPath.mkdir()
    get_NP_rings(inputPath, saveRingsPath)
    
    print('Finished in:')
    print(datetime.now() - start_time)
