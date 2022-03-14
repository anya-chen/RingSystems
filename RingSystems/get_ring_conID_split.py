import os
import pandas as pd
from rdkit import Chem
from pathlib import Path
import glob
from tqdm import tqdm
from datetime import datetime

'''
This script is used to give ring systems rign_conID and split files with the same ring_conID molecules in one file.(for further giving ring_stereoID.

use after get_NP_rings.py and get_ZINC_rings.py, get_organismSets_rings.py
'''

def assign_conID_split_files(inputfolder,outputfolder,prefix,size=150):
    ''' use for a big file has been splitted in parts.
    Need to put them together first and then split according to consititutions
    '''
    # get all .csv files in one folder
    csv_files = glob.glob(os.path.join(inputfolder, prefix+"*.csv"))
    frames = []
    for f in csv_files:
        frames.append(pd.read_csv(f,sep='\t'))
    df = pd.concat(frames)
    print('number of entries: {}'.format(str(len(df))))
    
    print('getting ringSmiles_noStereo......')
    df['ringSmiles_noStereo'] = df.apply(get_smiles_noStereo,axis=1)
    
    df_groupby = df.groupby('ringSmiles_noStereo',sort=False)
    frames = []
    conID = 0
    i = 0
    for name,group in df_groupby:
        group.count().reset_index()
        conID += 1
        group['ring_conID'] = conID
        frames.append(group)
        if len(frames) % size == 0:
            i += 1
            outputFileName = prefix + '_rings_p' + str(i) + '.csv'
            outputFile = outputfolder / outputFileName
            results_df = pd.concat(frames)
            del results_df['ringSmiles_noStereo']
            results_df.to_csv(outputFile,sep='\t',index=False)
            del results_df
            frames = []
            print('writing down part {} ......'.format(str(i)))
    i+=1
    outputFileName = prefix + '_rings_p' + str(i) + '.csv'
    outputFile = outputfolder / outputFileName
    results_df = pd.concat(frames)
    del results_df['ringSmiles_noStereo']
    results_df.to_csv(outputFile,sep='\t',index=False)
    del results_df
    print('writing down part {} ......'.format(str(i)))
    
    
def get_smiles_noStereo(row):
    smi = row['RingSmiles']
    smiles_noStereo = Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=False)
    return smiles_noStereo


if __name__ == "__main__":
    tqdm.pandas()
    start_time = datetime.now()
    
    refinedInputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/ringSystems')
    refinedOutputFolder = refinedInputFolder.parent / 'ringSystems_conID'
    if not refinedOutputFolder.exists():
        refinedOutputFolder.mkdir()

    assign_conID_split_files(refinedInputFolder,refinedOutputFolder,'refinedCoconut')

    
    zincInputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/ringSystems')
    zincOutputFolder = zincInputFolder.parent / 'ringSystems_conID'
    if not zincOutputFolder.exists():
        zincOutputFolder.mkdir()
        
    assign_conID_split_files(zincInputFolder,zincOutputFolder,'zinc')
    
    
    inputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/ringSystems')
    outputFolder = inputFolder.parent / 'ringSystems_conID'
    if not outputFolder.exists():
        outputFolder.mkdir()

    assign_conID_split_files(inputFolder,outputFolder,'plants',size=300)
    assign_conID_split_files(inputFolder,outputFolder,'bacteria',size=300)
    assign_conID_split_files(inputFolder,outputFolder,'fungi',size=300)
    assign_conID_split_files(inputFolder,outputFolder,'marine',size=300)
    
    print(datetime.now()-start_time)
    print('done')