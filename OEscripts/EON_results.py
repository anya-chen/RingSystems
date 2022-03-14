import os
import pandas as pd
from openeye import oechem
from datetime import datetime


'''
Script to obtain two different summarizing csv files of the EON results. One result file lists for each NP ring system
the best matching SM ring systems for each similarity measure individually. The others take only the best matching 
synthetic molecule ring system by ET_pb or ET_combo. 
Use after the application of EON.
'''


def save_best_eon_results(resultFilePath, eonDir, ZINCsmilesDict, NPsmilesDict, byScore):
    '''
    Saves the results of all EON result files in a tab separated file

    :param resultFilePath: the filepath to the tab separated result file
    :param eonDir: the directory to all EON result files
    :param ZINCsmilesDict: ID:SMILES dictionary of ZINC ring systems
    :param NPsmilesDict: ID:SMILES dictionary of NP ring systems
    :param byScore: for the results only the ring system best ranked by 'EON_ET_pb', 'EON_ET_combo' or 'all'
    '''
    resultDFlist = []
    for resultDict in get_eon_results(eonDir, ZINCsmilesDict, NPsmilesDict, byScore):
        resultDF = pd.DataFrame([resultDict])
        resultDFlist.append(resultDF)
    resultDF = pd.concat(resultDFlist, ignore_index=True)
    print(resultDF.info)

    resultDF.to_csv(resultFilePath, sep='\t', index=False, float_format='%.3f')


def get_eon_results(dir, ZINCsmilesDict, NPsmilesDict, byScore):
    '''
    Yields a result dictionary for each EON .rpt file in a directory. The result dictionary contains the SMILES and
    name of a query ring system and the SMILES, the name and the score of the best ranked synthetic molecule ring
    systems.

    :param dir: directory with EON .rpt files
    :param ZINCsmilesDict: name:SMIlES dictionary of the synthetic molecule ring systems
    :param NPsmilesDict: name:SMIlES dictionary of the natural product ring systems
    :param byScore: for the results only the ring system best ranked by 'EON_ET_pb', 'EON_ET_combo' or 'all'
    :return: the result dictionary for each file with the query and best ranked synthetic molecule ring systems
    '''
    fileList = get_file_list(dir)
    fileDict = dict()
    for file in fileList:
        file_prefix = '_'.join(file.split('_')[0:2])
        if file_prefix in fileDict.keys():
            fileDict[file_prefix].append(file)
        else:
            fileDict[file_prefix] = [file]
            
    resultList = ['EON_ET_pb', 'EON_ET_combo', 'EON_ShapeTanimoto']
    counter = 0
    for prefix,file_list in fileDict.items():
        frames = []
        for f in file_list:
            if counter % 1000 == 0:
                print('{} files processed.'.format(counter))
            
            filepath = os.path.join(dir, f)
            file_df = pd.read_csv(filepath, sep='\t')
            frames.append(file_df)
        eonResultDf = pd.concat(frames,ignore_index=True)

        if byScore == 'all':
            resultDict = get_best_results_for_resultList(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList)
        elif byScore =='EON_ET_pb':
            resultDict = get_best_eon_results_by_ET(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList)
        elif byScore == 'EON_ET_combo':
            resultDict = get_best_eon_results_by_combo(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList)
        counter += 1
        if resultDict:
            yield resultDict


def get_file_list(dir):
    '''
    Creates a file list with all files in a directory which end with .rpt

    :param dir: directory with files
    :return: a list of all .rpt files
    '''
    fileList = []
    for file in os.listdir(dir):
        if file.endswith('.rpt'):
            fileList.append(file)
    return fileList


def get_best_results_for_resultList(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList):
    '''
    gives the best results for each similarity measure in the given EON file

    :param eonResultDf: EON results for one natural product ring system stored in pandas dataframe
    :param NPsmilesDict: dictionary with smiles:id pairs of all natural product ring system SMILES
    :param ZINCsmilesDict: dictionary with smiles:id pairs of all ZINC ring system SMILES
    :param resultList: list with headers from the eon result file which should be evaluated
    :return: ring system IDs:SMILES dictionary for the ring systems with the best score for the given resultList entry
    '''
    resultDict = dict()
    resultDict['Query'] = eonResultDf.at[1, 'EON_Query']
    resultDict['Query_smiles'] = get_smiles_from_name(NPsmilesDict, resultDict['Query'])
    for resultListEntry in resultList:
        maxValue, maxIdx, ringSystemName = get_max_value(eonResultDf, resultListEntry)
        resultDict[resultListEntry] = maxValue
        resultDict[resultListEntry + "_name"] = ringSystemName
        resultDict[resultListEntry + "_smiles"] = get_smiles_from_name(ZINCsmilesDict, ringSystemName)
    if resultDict['EON_ET_pb'] == 1 and resultDict['Query_smiles'] == resultDict['EON_ET_pb_smiles']:
        return None
    return resultDict


def get_best_eon_results_by_ET(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList):
    '''
    gives the results for the best ranked synthetic ring system (ranked by EON_ET_pb)

    :param eonResultDf: EON results for one natural product ring system stored in pandas dataframe
    :param NPsmilesDict: dictionary with smiles:id pairs of all natural product ring system SMILES
    :param ZINCsmilesDict: dictionary with smiles:id pairs of all ZINC ring system SMILES
    :return:ring system IDs:SMILES dictionary for the best ranked synthetic ring system
    '''
    resultDict = dict()
    resultDict['Query'] = eonResultDf.at[1, 'EON_Query']
    resultDict['Query_smiles'] = get_smiles_from_name(NPsmilesDict, resultDict['Query'])
    maxValue, maxIdx, ringsystemName = get_max_value(eonResultDf, 'EON_ET_pb')
    resultDict['Best_match_by_ET'] = ringsystemName
    resultDict['Best_match_by_ET_SMILES'] = get_smiles_from_name(ZINCsmilesDict, ringsystemName)
    eonDict = eonResultDf.to_dict(orient='index')
    maxEonEntryDict = eonDict[maxIdx]
    for resultListEntry in resultList:
        resultDict[resultListEntry] = maxEonEntryDict[resultListEntry]

    if resultDict['EON_ET_pb'] == 1 and resultDict['Query_smiles'] == resultDict['Best_match_by_ET_SMILES']:
        return None
    return resultDict


def get_best_eon_results_by_combo(eonResultDf, NPsmilesDict, ZINCsmilesDict, resultList):
    '''
    gives the results for the best ranked synthetic ring system (ranked by EON_ET_combo)

    :param eonResultDf: EON results for one natural product ring system stored in pandas dataframe
    :param NPsmilesDict: dictionary with smiles:id pairs of all natural product ring system SMILES
    :param ZINCsmilesDict: dictionary with smiles:id pairs of all ZINC ring system SMILES
    :return:ring system IDs:SMILES dictionary for the best ranked synthetic ring system
    '''
    resultDict = dict()
    resultDict['Query'] = eonResultDf.at[1, 'EON_Query']
    resultDict['Query_smiles'] = get_smiles_from_name(NPsmilesDict, resultDict['Query'])
    maxValue, maxIdx, ringsystemName = get_max_value(eonResultDf, 'EON_ET_combo')
    resultDict['Best_match_by_combo'] = ringsystemName
    resultDict['Best_match_by_combo_SMILES'] = get_smiles_from_name(ZINCsmilesDict, ringsystemName)
    eonDict = eonResultDf.to_dict(orient='index')
    maxEonEntryDict = eonDict[maxIdx]
    for resultListEntry in resultList:
        resultDict[resultListEntry] = maxEonEntryDict[resultListEntry]

    if resultDict['EON_ET_combo'] == 1 and resultDict['Query_smiles'] == resultDict['Best_match_by_combo_SMILES']:
        return None
    return resultDict


def get_max_value(dataframe, columnName):
    '''
    Returns the maximum value of a given column in a dataframe of an EON report file

    :param dataframe: dataframe of an EON report file
    :param columnName: column for which the maximum value is requested
    :return: maximum value, the index and the name of the ring system with the maximum value
    '''
    column = dataframe[columnName]
    maxValue = column.max()
    maxIdx = column.idxmax()
    ringsystemName = dataframe.at[maxIdx, 'Name']
    name_split = ringsystemName.split('_')
    if len(name_split) == 3:
        ringsystemName = '_'.join(name_split[0:2])
    return maxValue, maxIdx, ringsystemName


def get_smiles_from_name(smilesDict, name):
    return smilesDict[name]


def get_smilesDic(file):
    ringDF = pd.read_csv(file,header=None,sep='\t')
    smilesDict = dict()
    for a,b in zip(ringDF[1],ringDF[0]):
        if a in smilesDict.keys():
            smilesDict[a] = smilesDict[a]+'.'+b
        else:
            smilesDict[a] = b
    return smilesDict
    


if __name__ == "__main__":
    start_time = datetime.now()
    
    npRings = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueRingSystems/np_uniqueRingSystems_for_omega.smi'
    zincRings = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/uniqueRingSystems/zinc_uniqueRingSystems_for_omega.smi'
    
    NPsmilesDict = get_smilesDic(npRings)
    ZINCsmilesDict = get_smilesDic(zincRings)

    save_best_eon_results('/data/local/ringsys/202111_3dComparison/results/EON_result_summary_by_best_ET.csv',
                          '/data/local/ringsys/202111_3dComparison/eon', ZINCsmilesDict, NPsmilesDict, 'EON_ET_pb')
    
    save_best_eon_results('/data/local/ringsys/202111_3dComparison/results/EON_result_summary_by_best_combo.csv',
                          '/data/local/ringsys/202111_3dComparison/eon', ZINCsmilesDict, NPsmilesDict, 'EON_ET_combo')
    
    save_best_eon_results('/data/local/ringsys/202111_3dComparison/results/EON_result_summary.csv',
                          '/data/local/ringsys/202111_3dComparison/eon', ZINCsmilesDict, NPsmilesDict,'all')
    
    print("Finished in:")
    print(datetime.now()-start_time)