import os
import multiprocessing as mp
from rdkit import Chem 

'''
Script to remove the NPs from the SMs by using the tautomer standardized SMILES without stereochemistry (after removing compounds with bad patterns, remove duplicates). It also creates unique NP files by preprocessed SMILES

Use after remove_duplicate*.py 
'''

def get_unique_set(NpDir):
    '''
    Create a set of unique tautomer standardized SMILES without stereochemistry
    
    :param NpDir: Input directory with preprocessed and tautomer standardized SMILES. The files are tab separated and have
    the columns preprocessedSmiles, id, conID, tautomerizedSmiles, with preprocessedSmiles being the SMILES without
    the tautomer standardization step.
    :return: update a set with tautomer standardized SMILES (without stereochemistry) to remove NP molecules from the ZINC set.
    '''    
    for file in os.listdir(NpDir):
        inputFilePath = os.path.join(NpDir, file)
        with open(inputFilePath, 'r', encoding='utf-8') as smilesCSV:
            header = smilesCSV.readline()
            for line in smilesCSV:
                entry = line.rstrip().split('\t')
                if len(entry) == 4:
                    tautomerizedSmiles = entry[3]
                    smiles_noStereo = Chem.MolToSmiles(Chem.MolFromSmiles(tautomerizedSmiles),isomericSmiles=False)
                    npTautoSmilesSetForZINC.add(smiles_noStereo)
        smilesCSV.close()
        print("Finished processing {}.".format(file))
    return npTautoSmilesSetForZINC


def remove_NP_from_preprocessed_ZINC_smiles(npTautoSmilesSetForZINCWithCoconut, processedZINCSmilesPath, outputFilepath):
    '''
    Remove NPs in the whole coconut for ZINC
    
    :param npTautoSmilesSetForZINCWithCoconut: NP SMILES set with all tautomer standardized SMILES including removedCoconut molecules
    :param processedZINCSmilesPath: filepath to ZINC file with preprocessed ZINC molecules as SMILES. The files are tab
    separated and have the columns preprocessedSmiles, id, conID, tautomerizedSmiles.
    :param outputFilepath: filepath where the ZINC molecules are saved after removing NPs
    '''
    with open(processedZINCSmilesPath, 'r') as ZINCsmilesCSV:
        with open(outputFilepath, 'w') as cleanedZINCfile:
            header = ZINCsmilesCSV.readline()
            cleanedZINCfile.write(header)
            for line in ZINCsmilesCSV:
                entry = line.rstrip().split('\t')
                if len(entry) == 4:
                    tautomerizedSmiles = entry[3]
                    smiles_noStereo = Chem.MolToSmiles(Chem.MolFromSmiles(tautomerizedSmiles),isomericSmiles=False)
                    if smiles_noStereo not in npTautoSmilesSetForZINCWithCoconut:
                        cleanedZINCfile.write(line)
        cleanedZINCfile.close()
    ZINCsmilesCSV.close()
    print('{} Successfully cleaned from NPs.'.format(processedZINCSmilesPath))


if __name__ == "__main__":

    refinedCoconutDir = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueSmiles'
    removedCoconutDir = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut/uniqueSmiles'
    zincInputDir = '/home/ychen/projects/nps_ringsys/20210519_data_prep/ZINC/uniqueSmiles'
    zincOutputDir = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/uniqueSmilesCleaned'

    if not os.path.exists(zincOutputDir):
        os.makedirs(zincOutputDir)

    npTautoSmilesSetForZINC = set()
    npTautoSmilesSetForZINC = get_unique_set(refinedCoconutDir)
    print("NP Set completed holding {} molecules.".format(len(npTautoSmilesSetForZINC)))

    npTautoSmilesSetForZINCWithRemovedCoconut = get_unique_set(removedCoconutDir)
    print("removedCoconut added to NP Set, holding {} molecules.".format(len(npTautoSmilesSetForZINCWithRemovedCoconut)))

    jobList = []
    for file in os.listdir(zincInputDir):
        processedZINCSmilesPath = os.path.join(zincInputDir, file)
        outputFile = "NP_removed_" +file
        outputSmilesFilePath = os.path.join(zincOutputDir, outputFile)
        jobList.append((npTautoSmilesSetForZINCWithRemovedCoconut, processedZINCSmilesPath, outputSmilesFilePath))

    pool = mp.Pool(processes=mp.cpu_count()-1)
    pool.starmap(remove_NP_from_preprocessed_ZINC_smiles, jobList)
    pool.close()
    pool.join()
