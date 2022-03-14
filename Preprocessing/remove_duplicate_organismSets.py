import glob
from pathlib import Path
from rdkit import Chem
import pandas as pd
from datetime import datetime


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



if __name__ == "__main__":
    start_time = datetime.now()
    
    inputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/noBadPatternSmiles')
    outputFolder = inputFolder.parent / 'uniqueSmiles'
    if not outputFolder.exists():
        outputFolder.mkdir()
    
    files = inputFolder.glob('*')
    for file in files:
        print(file)
        df = pd.read_csv(file)
        df = df[['preprocessedSmiles','id','tautomerizedSmiles']]
        df = df[df['preprocessedSmiles'] != '']
        print('number of entries after remove empty preprocessedSmiles: {}'.format(str(len(df))))
        
        df.drop_duplicates(subset=['preprocessedSmiles'],inplace=True)
        print('number of entries after drop duplicated preprocessedSmiles: {}'.format(str(len(df))))
    
        df['smiles_noStereo'] = df.apply(get_smiles_noStereo,axis=1)
        df['nChiralFlags'] = df.apply(count_nFlags_smiles,axis=1)
        
        df = df.groupby('smiles_noStereo')
        
        frames = []
        conID = 0
        i = 0
        for name,group in df:
            group.count().reset_index()
            conID += 1
            group['conID'] = conID  # can not compare this ID of compounds in different datasets.
            
            nums = set(group.nChiralFlags)
            if len(nums) >1 and (0 in nums):
                new_group = group[group.nChiralFlags != 0]
            else:
                new_group = group

            new_group.count().reset_index()
            frames.append(new_group)
            if len(frames) % 3000 == 0:
                i +=1
                results_df = pd.concat(frames)
                outFileName = str(file.stem)+'_uniq_'+str(i)+'.csv'
                outputFile = outputFolder / outFileName
                results_df[['preprocessedSmiles','id','conID','tautomerizedSmiles']].to_csv(outputFile,sep='\t',index=False)
                frames = []
        i +=1
        results_df = pd.concat(frames)
        outFileName = str(file.stem)+'_uniq_'+str(i)+'.csv'
        outputFile = outputFolder / outFileName
        results_df[['preprocessedSmiles','id','conID','tautomerizedSmiles']].to_csv(outputFile,sep='\t',index=False)
    
    print(datetime.now()-start_time)
    print('Done')