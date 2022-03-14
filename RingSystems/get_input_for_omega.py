import os
import glob
from pathlib import Path
import pandas as pd
from datetime import datetime

'''
This script is used to get the unique smiles of ring systems as input for omega.
output as .smi with smiles id

run after RingSystems/assign_ring_stereoIDs.py
'''
    

def get_unique_ringsys_for_omega(inputFolder,outFile):
    csv_files = glob.glob(os.path.join(inputFolder, "*.csv"))
    dframes = []
    for f in csv_files:
        df = pd.read_csv(f)
        df_stereo = unique_ringsys_for_omega(df)
        dframes.append(df_stereo)
    results_df = pd.concat(dframes) 
    results_df.to_csv(outFile,sep='\t',columns=['RingSmiles','ring_stereoIDs'],header=False, index=False)

    
def unique_ringsys_for_omega(df):
    df.drop_duplicates('RingSmiles',keep='first',inplace=True)
    
    df_stereoID = df[['RingSmiles','ring_stereoIDs']]
    df_stereoID['nFlags'] = df_stereoID.apply(lambda row:count_nFlags_smiles(row.RingSmiles),axis=1)
    df_stereoID['nbFlags'] = df_stereoID.apply(lambda row:count_bondFlags_smiles(row.RingSmiles),axis=1)
                    
    # remove rows with list of ids
    frames = []
    ids = []
    for row in df_stereoID.itertuples():
        if '[' not in row.ring_stereoIDs:
            frames.append(row)

    df_stereoID_sep = pd.DataFrame(frames)
    del df_stereoID_sep['Index']
    
    df_stereoID_groups = df_stereoID_sep.groupby('ring_stereoIDs')
    frames = []
    for name,group in df_stereoID_groups:
        if len(group)==1:
            frames.append(group)
        else:
            new_group = group[group['nFlags'] == sorted(set(group.nFlags))[-1]]
            new_group = new_group[new_group['nbFlags'] == sorted(set(group.nbFlags))[-1]]
            if len(new_group) != 1 and set(new_group['nbFlags']) == {0}:
                print(new_group)
                new_group.drop_duplicates('ring_stereoIDs',inplace=True)
            frames.append(new_group)
    df_stereo = pd.concat(frames)
    return df_stereo


def count_nFlags_smiles(smi):
    smi = smi.replace('@@','@')
    nFlags = smi.count('@')
    return nFlags


def count_bondFlags_smiles(smi):
    smi = smi.replace('\\','/')
    nbFlags = smi.count('/')
    return nbFlags

    
    
if __name__ == "__main__":
    start_time = datetime.now()
    
    #np
    npFolder = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/ringSystems_stereoID'
    npOutputDir  = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueRingSystems')
    if not npOutputDir.exists():
        npOutputDir.mkdir()
    npOutputFile = npOutputDir / 'np_uniqueRingSystems_for_omega.smi'
    
    get_unique_ringsys_for_omega(npFolder,npOutputFile)
    
    #zinc
    zincFolder = npFolder.replace('refinedCoconut','zinc')
    zincOutputDir = Path(str(npOutputDir).replace('refinedCoconut','zinc'))
    if not zincOutputDir.exists():
        zincOutputDir.mkdir()
    zincOutputFile = zincOutputDir / 'zinc_uniqueRingSystems_for_omega.smi'
    
    get_unique_ringsys_for_omega(zincFolder,zincOutputFile)
    print(datetime.now()-start_time)
    print('done')
