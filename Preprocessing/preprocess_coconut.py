import os
from pathlib import Path
import csv
from rdkit import Chem
from split_NP_files import split_large_files
from preprocess_molecules import preprocess_database
import multiprocessing as mp

'''
The script to process COCONUT DB, including common preprocessing and tautomerization in order to remove potential natural products from the ZINC dataset. 
Use before remove_NPs_from_ZINC.py
'''

def get_raw_data_from_coconut(inputFile,reformatedFile):
    '''
    Reads in the file of the COCONUT DB and returns a ID:SMILES dictionary.

    :param inputfile: File with SMILES and IDs obtained from COCONUT
    :return: a ID:SMILES dictionary
    '''
    
    smilesDict = dict()
    with open(reformatedFile,'w',encoding='utf-8') as f:
        with open(inputFile, 'r', encoding='utf-8') as coconutFile:
            moleculeCsv = csv.reader(coconutFile, delimiter='\t')
            next(moleculeCsv)
            for row in moleculeCsv:
                if len(row)<8:continue
                smiles = row[2]
                # the following lines were used to remove atom numbers in the input smiles
                # also remove smiles cannot convert to molecules with rdkit
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
                    smiles = Chem.MolToSmiles(mol)
                    name = row[0]
                    f.write('{}\t\t= {}\n'.format(name,smiles))
                    smilesDict[name] = smiles
            
    return smilesDict


def get_molecule_dict_as_smi_output(smilesDict, outputFile):
    '''
    Writes the ID:SMILES dict in a new file.

    :param smilesDict: ID:SMILES dictionary
    :param outputFile: filepath to smiles file
    '''
    with open(outputFile, 'w', encoding="utf-8") as of:
        of.write('smiles ID\n')
        for moleculeName, smiles in smilesDict.items():
            if smiles != '':
                of.write('{} {}\n'.format(smiles, moleculeName))
    
def preprocess(inputFile,outputFolder):
    if not outputFolder.exists():
        outputFolder.mkdir()
    
    #Reformat COCONUT_DB.smi
    reformatedFolder = outputFolder / 'smiles_reformat'
    if not reformatedFolder.exists():
        reformatedFolder.mkdir()    
    reformatedName = inputFile.split('/')[-1].split('.')[0]+'.smi'
    reformatedFile = reformatedFolder / reformatedName
                         
    smilesDict = get_raw_data_from_coconut(inputFile,reformatedFile)
    get_molecule_dict_as_smi_output(smilesDict, reformatedFile)
    
    #Split it into multiple files
    splitFolder = outputFolder / 'smiles_input'
    if not splitFolder.exists():
        splitFolder.mkdir() 
    
    split_large_files(reformatedFile,3000)
    
    #structure preprocess
    preprocessedDir = outputFolder / 'preprocessedSmiles'
    logdir = outputFolder / 'log'
    if not preprocessedDir.exists():
        preprocessedDir.mkdir()
    if not logdir.exists():
        logdir.mkdir()
    
    databaseFilesJobList = []
    for file in os.listdir(splitFolder):
        databaseFilesJobList.append(os.path.join(splitFolder, file))
    pool = mp.Pool(processes=mp.cpu_count()-1)
    pool.map(preprocess_database, databaseFilesJobList)
    pool.close()
    pool.join()
                    
if __name__ == "__main__":

    refinedInputFile = "/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/refinedCoconut.sourceNP.csv"
    refinedOutputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut')
    preprocess(refinedInputFile,refinedOutputFolder)
    
    removedInputFile = "/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/removedCoconut.sourceNP.csv"
    removedOutputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut')
    preprocess(removedInputFile,removedOutputFolder)
    
    print('Done')
    
