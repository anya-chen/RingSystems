import os
from pathlib import Path
from src import MoleculePreprocessorExtended
from rdkit import RDLogger
from rdkit import Chem

'''
Preprocess molecules using the MoleculePreprocessor Class, expects SMILES files as input
Use for all molecules, NPs and SMs as well.
Used in reprocess_Zinc.py and preprocess_coconut.py
'''


def preprocess_database(file):
    '''
    processes the molecules of a file and save them as a csv in ../preprocessedSmiles/outFilename.csv
    :param file: file with molecules to be processed.
    '''
    fileName = os.path.basename(file).split('.')[0]
    outFilename = fileName + '_preprocessed.csv'
    outpath = Path(file).parent.parent / 'preprocessedSmiles' / outFilename
        
    if Path(outpath).is_file():
        raise Exception("File {} already exists!".format(outpath))

    smilesIDdict = convert_csv_to_smiles_dict(file)
    print('{} successfully converted to smilesDict'.format(file))

    listOfSmiles = list(smilesIDdict.keys())

    rdLogger = RDLogger.logger()

    moleculesProcessed = MoleculePreprocessorExtended.MoleculePreprocessorExtended.init_with_smiles(listOfSmiles)
    
    # set c++ log level of rdkit temporarily to ERROR so that csp does not clutter
    # stdout with logging
    rdLogger.setLevel(RDLogger.ERROR)
    moleculesProcessed.csp()  #use the chembl_structure_pipeline (standardizer and get_parent)
    rdLogger.setLevel(RDLogger.INFO)
    
    moleculesProcessed.element_filter()
    moleculesProcessed.check_molecules_validity()
    # rawsmiles:preprocessedSmiles dict
    preprocessedSmilesDict = moleculesProcessed.get_rawsmiles_smiles_dict()
    
    print('{} successfully preprocessed.'.format(file))
    
    moleculesProcessed.canonalize_tautomer()
    # rawsmiles:tautomerizedSmiles dict
    tautoMoleculesDict = moleculesProcessed.get_rawsmiles_smiles_dict()

    print('{} successfully tautomerized.'.format(file))

    logFileName = fileName + '.log'
    logFilepath = Path(file).parent.parent / 'log' / logFileName
    moleculesProcessed.write_log_file(logFilepath)

    #(rawSmiles, preprocessedSmiles, tautomerizedSmiles, ID)
    finalDict = mergeDict(tautoMoleculesDict, mergeDict(smilesIDdict, preprocessedSmilesDict))
    
    save_dict_as_csv(finalDict, outpath, ['rawSmiles', 'preprocessedSmiles', 'id', 'tautomerizedSmiles'])


def convert_csv_to_smiles_dict(file):
    '''
    expects a csv file with header line (which is skipped) and returns a rawsmiles:databaseID dict

    :param file: a csv file, using 'space' as delimiter
    :return: a rawsmiles:databaseID dict
    '''
    # rawsmiles:databaseID dict
    smilesIDdict = dict()
    with open(file, 'r', encoding='utf8') as smilesIDfile:
        # skip first (header) line
        next(smilesIDfile)
        for line in smilesIDfile:
            entry = line.rstrip().split()
            smiles = entry[0]
            databaseID = entry[1]
            smilesIDdict[smiles] = databaseID
        return smilesIDdict


def mergeDict(dict1, dict2):
    ''' 
    Merge dictionaries and keep values of common keys in list
    '''
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = [value, dict1[key]]
    return dict3


def save_dict_as_csv(mergedDict, outputFile, headerList):
    '''
    Saves the preprocessed molecules as smiles in a tab seperated csv

    :param mergedDict: a rawSmiles: [[prepSmiles, id], canonSmiles] dictionary
    :param outputFile: filepath to the output file
    :param headerList: list with headers for the csv
    '''
    
    with open(outputFile, 'w', encoding='utf8') as csvFile:
        csvFile.write('\t'.join(headerList) + '\n')
        for rawsmiles, values in mergedDict.items():
            csvFile.write('{}\t{}\t{}\t{}\n'.format(rawsmiles, values[0][0], values[0][1], values[1]))
    print('{} successfully saved'.format(outputFile))

