import os
import glob
import pandas as pd
from pathlib import Path
from rdkit import Chem


'''
This script is used to get the unique smiles of ring systems, both with and without stereochemistry.
also count how many compounds have a certain ring system.
the results has also full list of ids of those compounds and ring systmes incase nned to look up.
+ save two other simplier files: 
    .txt with smiles,id,nMol
    .smi with smiles,id

run after RingSystems/assign_ring_stereoIDs.py
          RingSystems/get_organismSets_ids.py
'''

def get_unique_ringsys_noStereo(inputFolder,outFile_noStereo,prefix):
    csv_files = glob.glob(os.path.join(inputFolder, prefix+"*.csv"))
    frames = []
    for f in csv_files:
        df = pd.read_csv(f)
        df_noStereo = unique_ringsys_noStereo(df)
        frames.append(df_noStereo)
        
    results_df = pd.concat(frames)
    save_unique_ringsys_noStereo(results_df,outFile_noStereo)
        

def unique_ringsys_noStereo(df):
    df_conID = df[['RingSmiles','ring_conID','MoleculeId','MoleculeConID']] 
    groupped_df_conID = df_conID.groupby('ring_conID')
    df_noStereo = groupped_df_conID.agg({'MoleculeConID':'unique','MoleculeId':'unique',
                                         'RingSmiles':'first'}) #'MoleculeConID':'nunique'
    df_noStereo = df_noStereo.reset_index()
    df_noStereo['nMol_conID'] = df_noStereo.apply(lambda row:len(row.MoleculeConID),axis=1)
    df_noStereo['ringSmiles_noStereo'] = df_noStereo.apply(lambda row:Chem.MolToSmiles(Chem.MolFromSmiles(row.RingSmiles),isomericSmiles=False),axis=1)
    del df_noStereo['RingSmiles']
    return df_noStereo


def save_unique_ringsys_noStereo(results_df,outFile_noStereo):
    results_df.sort_values('nMol_conID',ascending=False,inplace=True)
    results_df.to_csv(outFile_noStereo,sep='\t',columns=['ringSmiles_noStereo','ring_conID','nMol_conID','MoleculeConID','MoleculeId'],index=False)
    outFile_txt = str(outFile_noStereo).split('.')[0] + '.txt'
    results_df.to_csv(outFile_txt,sep='\t',columns=['ringSmiles_noStereo','ring_conID','nMol_conID'],index=False)
    outFile_smi = str(outFile_noStereo).split('.')[0] + '.smi'
    results_df.to_csv(outFile_smi,sep='\t',columns=['ringSmiles_noStereo','ring_conID'],index=False)
    

def get_unique_ringsys_Stereo(inputFolder,outFile_Stereo,prefix):        
    csv_files = glob.glob(os.path.join(inputFolder, prefix+"*.csv"))
    dataframes = []
    for f in csv_files:
        df = pd.read_csv(f)
        df_Stereo = unique_ringsys_Stereo(df)
        dataframes.append(df_Stereo)
        
    results_df = pd.concat(dataframes)
    save_unique_ringsys_Stereo(results_df,outFile_Stereo)
    
    
def unique_ringsys_Stereo(df):
    df_stereoID = df[['RingSmiles','ring_conID','ring_stereoIDs','MoleculeId','MoleculeConID','mol_stereoIDs']]
    df_stereoID['nFlags'] = df_stereoID.apply(lambda row:count_nFlags_smiles(row.RingSmiles),axis=1)
      
    # separate multi-ids rings as multiple rows
    frames = []
    ids = []
    for row in df_stereoID.itertuples():
        if '[' in row.ring_stereoIDs:
            id_str = row.ring_stereoIDs.replace('[','').replace(']','').replace("'","").replace(' ','')
            id_list = id_str.split(',')
            for i in id_list:
                ids.append(i)
                frames.append(row)
        else:
            ids.append(row.ring_stereoIDs)
            frames.append(row)

    df_stereoID_sep = pd.DataFrame(frames)
    del df_stereoID_sep['Index']
    df_stereoID_sep['ring_stereoID_sep'] = ids
    
    # drop rows with ‘nFlags’ not equal to the largest available. (feel can just drop all multi-ids rings in the last step, but still need to separate ids with '[' and do this step anyways)
    groupped = df_stereoID_sep.groupby('ring_stereoID_sep')
    frames = []
    for name,group in groupped:
        if len(group) == 1:
            frames.append(group)
        else:
            new_group = group[group['nFlags'] == sorted(set(group.nFlags))[-1]]
            frames.append(new_group)
    df_stereoID_sep = pd.concat(frames)

    # keep the first smiles of each ring_stereoID
    df_stereoID_unique = df_stereoID_sep[['RingSmiles','ring_conID','ring_stereoID_sep']]
    df_stereoID_unique.drop_duplicates('ring_stereoID_sep',keep='first',inplace=True)

    # separate multi-ids molecules as multiple rows
    frames = []
    ids = []
    for row in df_stereoID_sep.itertuples():
        if '[' in row.mol_stereoIDs:
            id_str = row.mol_stereoIDs.replace('[','').replace(']','').replace("'","").replace(' ','')
            id_list = id_str.split(',')
            for i in id_list:
                ids.append(i)
                frames.append(row)
        else:
            ids.append(row.mol_stereoIDs)
            frames.append(row)

    df_stereoID_sep2 = pd.DataFrame(frames)
    del df_stereoID_sep2['Index']
    df_stereoID_sep2['mol_stereoID_sep'] = ids

    # now can count
    groupped_df_stereoID = df_stereoID_sep2.groupby('ring_stereoID_sep')
    df_stereoID_count = groupped_df_stereoID.agg({'MoleculeConID':'unique',
                                                  'MoleculeId':'unique',
                                                 'mol_stereoID_sep':'unique'}) 
    df_stereoID_count = df_stereoID_count.reset_index()
    df_stereoID_count['nMol_stereoID'] = df_stereoID_count.apply(lambda row:len(row.mol_stereoID_sep),axis=1)
    df_stereo = pd.merge(df_stereoID_unique,df_stereoID_count,on='ring_stereoID_sep')
    df_stereo.rename(columns={'ring_stereoID_sep':'ring_stereoID',
                      'mol_stereoID_sep':'mol_stereoID_list'},inplace = True)
    return df_stereo
  
    
def count_nFlags_smiles(smi):
    smi = smi.replace('@@','@')
    nFlags = smi.count('@')
    return nFlags

    
def save_unique_ringsys_Stereo(results_df,outFile_Stereo):
    results_df.sort_values('nMol_stereoID',ascending=False,inplace=True)
    results_df.to_csv(outFile_Stereo,sep='\t',columns=['RingSmiles','ring_stereoID','nMol_stereoID','ring_conID','MoleculeId','MoleculeConID','mol_stereoID_list'],index=False)
    outFile_txt = str(outFile_Stereo).split('.')[0] + '.txt'
    results_df.to_csv(outFile_txt,sep='\t',columns=['RingSmiles','ring_stereoID','nMol_stereoID'],index=False)
    outFile_smi = str(outFile_Stereo).split('.')[0] + '.smi'
    results_df.to_csv(outFile_smi,sep='\t',columns=['RingSmiles','ring_stereoID'],index=False)

    
    
if __name__ == "__main__":
    # refinedCoconut
    npInputDir = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/ringSystems_stereoID'
    npOutputDir  = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueRingSystems')
    if not npOutputDir.exists():
        npOutputDir.mkdir()
    npOutput1 = npOutputDir / 'np_uniqueRingSystems_noStereo.csv'
    npOutput2 = npOutputDir / 'np_uniqueRingSystems_Stereo.csv'
                     
    get_unique_ringsys_noStereo(npInputDir,npOutput1,'refinedCoconut')
    get_unique_ringsys_Stereo(npInputDir,npOutput2,'refinedCoconut')
    
    # zinc
    zincInputDir = npInputDir.replace('refinedCoconut','zinc')
    zincOutputDir = Path(str(npOutputDir).replace('refinedCoconut','zinc'))
    if not zincOutputDir.exists():
        zincOutputDir.mkdir()
    zincOutput1 = zincOutputDir / 'zinc_uniqueRingSystems_noStereo.csv'
    zincOutput2 = zincOutputDir / 'zinc_uniqueRingSystems_Stereo.csv'
                     
    get_unique_ringsys_noStereo(zincInputDir,zincOutput1,'zinc')
    get_unique_ringsys_Stereo(zincInputDir,zincOutput2,'zinc')
    
    # orgainismSets
    subsetsDir = npInputDir.replace('refinedCoconut','organismSets')
    subsets_OutputDir = Path(str(npOutputDir).replace('refinedCoconut','organismSets'))
    if not subsets_OutputDir.exists():
        subsets_OutputDir.mkdir()
    
    subset_list = ['plants','bacteria','fungi','marine']
    for subset in subset_list:
        outName1 = subset+'_uniqueRingSystems_noStereo.csv'
        outName2 = subset+'_uniqueRingSystems_Stereo.csv'
        subsetsOutput1 = subsets_OutputDir / outName1
        subsetsOutput2 = subsets_OutputDir / outName2
        get_unique_ringsys_noStereo(subsetsDir,subsetsOutput1,subset)
        get_unique_ringsys_Stereo(subsetsDir,subsetsOutput2,subset)
    
    
    print('done')
