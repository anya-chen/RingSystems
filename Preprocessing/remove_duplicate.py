import os
import pandas as pd
from rdkit import Chem
from pathlib import Path


'''
This script is used to remove duplicate molecules for coconut, zinc (datasets have been splitted in parts).
It will:
    - group by unique constitution (same connectivity) and have a constitution ID (conID).
    - within one group, remove smiles with fully undefined stereochemistry if there is a better defined smiles (only consider @ flags in smiles)
    
use after preprocessing, filter bad patterns.

'''


def remove_duplicate_splitted_files(inputfolder,outputfolder,prefix):
    ''' use for a big file has been splitted in parts.
    Need to put them together first and then split according to consititutions
    '''
    li = []
    for filename in os.listdir(inputfolder):
        file = inputfolder / filename
        df = pd.read_csv(file,sep='\t',index_col=None, header=0)
        li.append(df)

    df = pd.concat(li, axis=0, ignore_index=True)
    print('number of entries: {}'.format(str(len(df))))

    df = df[['preprocessedSmiles','id','tautomerizedSmiles']]
    df = df[df['preprocessedSmiles'] != '']
    print('number of entries after remove empty preprocessedSmiles: {}'.format(str(len(df))))

    df.drop_duplicates(subset=['preprocessedSmiles'],inplace=True)
    print('number of entries after drop duplicated preprocessedSmiles: {}'.format(str(len(df))))
    
    df['smiles_noStereo'] = df.apply(get_smiles_noStereo,axis=1)
    df['nChiralFlags'] = df.apply(count_nFlags_smiles,axis=1)
    
    df_noVariance,df_Variance = separate_by_if_stereo_variance(df)
    del df
    
    outputfilename1 = prefix + '_noVariance_uniq.csv'
    outputFile1 = outputfolder/ outputfilename1
    df_noVariance[['preprocessedSmiles','id','conID','tautomerizedSmiles']].to_csv(outputFile1,sep='\t',index=False)
    print('number of entries with no stereo Variance: {}'.format(str(len(df_noVariance))))
    del df_noVariance
    
    variance = len(df_Variance)
    print('number of entries with stereo Variance: {}'.format(str(variance)))
    
    outputfilename2 = prefix + '_Variance_uniq.csv'
    outputFile2 = outputfolder/ outputfilename2
    split_dataframe_by_groups_with_tautomerizedSmiles(df_Variance,outputFile2)

    
def get_smiles_noStereo(row):
    smi = row['preprocessedSmiles']
    smiles_noStereo = Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=False)
    return smiles_noStereo


def count_nFlags_smiles(row):
    #annotate how many chiral tags in the smiles
    smi = row['preprocessedSmiles']
    smi = smi.replace('@@','@')
    nFlags = smi.count('@')
    return nFlags


def separate_by_if_stereo_variance(df):
    df = df.groupby('smiles_noStereo')
    
    frames = []
    frames_variance = []
    conID = 0
    for name,group in df:
        group.count().reset_index()
        conID += 1
        group['conID'] = conID #can not compare this ID of compounds in different datasets.
        
        if len(group) == 1:
            frames.append(group)
        else:
            frames_variance.append(group)
            
    df_noVariance = pd.concat(frames)
    df_Variance = pd.concat(frames_variance)
    return df_noVariance,df_Variance

    
def split_dataframe_by_groups_with_tautomerizedSmiles(df,outputFile,size=5000):
    '''use this to split a set to smaller parts
    '''
    df = df.groupby('smiles_noStereo')

    frames = []
    i = 0
    for name,group in df:
        group.count().reset_index()
        frames.append(group)
        if len(frames) % size == 0:
            i += 1
            results_df = remove_smiles_fully_undefined_stereo(pd.concat(frames))
            outputFile3 = str(outputFile).split('.')[0] + '_part'+str(i)+'.csv' 
            print('writing down part {} ......'.format(str(i)))
            results_df[['preprocessedSmiles','id','conID','tautomerizedSmiles']].to_csv(outputFile3,sep='\t',index=False)
            del results_df
            frames = []
    i += 1
    results_df = remove_smiles_fully_undefined_stereo(pd.concat(frames))
    outputFile3 = str(outputFile).split('.')[0] + '_part'+str(i)+'.csv' 
    print('writing down part {} ......'.format(str(i)))
    results_df[['preprocessedSmiles','id','conID','tautomerizedSmiles']].to_csv(outputFile3,sep='\t',index=False)
    del results_df


def remove_smiles_fully_undefined_stereo(df): 
    '''only when there is other smiles in this consititution
    '''
    df = df.groupby('smiles_noStereo')    
    
    frames = []
    for name,group in df:
        group.count().reset_index()
        nums = set(group.nChiralFlags)
        
        if len(nums) >1 and (0 in nums):
            new_group = group[group.nChiralFlags != 0]
        else:
            new_group = group
                
        new_group.count().reset_index()
        frames.append(new_group)

    results_df = pd.concat(frames)
    return results_df


    
if __name__ == "__main__":
    refinedInputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/noBadPatternSmiles')
    refinedOutputFolder = refinedInputFolder.parent / 'uniqueSmiles'
    if not refinedOutputFolder.exists():
        refinedOutputFolder.mkdir()

    remove_duplicate_splitted_files(refinedInputFolder,refinedOutputFolder,'refinedCoconut')


    removedInputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut/noBadPatternSmiles')
    removedOutputFolder = removedInputFolder.parent / 'uniqueSmiles'
    if not removedOutputFolder.exists():
        removedOutputFolder.mkdir()

    remove_duplicate_splitted_files(removedInputFolder,removedOutputFolder,'removedCoconut')
    
    
    zincInputFolder = Path('/home/ychen/projects/nps_ringsys/20210519_data_prep/ZINC/noBadPatternSmiles')
    zincOutputFolder = zincInputFolder.parent / 'uniqueSmiles'
    if not zincOutputFolder.exists():
        zincOutputFolder.mkdir()
        
    remove_duplicate_splitted_files(zincInputFolder,zincOutputFolder,'zinc')
    
    print('Done')
