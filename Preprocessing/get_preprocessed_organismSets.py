import os
from pathlib import Path
import glob
import pandas as pd
from rdkit import Chem

'''
This script is used to get preprocessed NP subsets from different organisms
since the whole refinedCoconut is preprocessed, here should only use the counter list to get directly from the preprocessed file of refinedCoconut. 
#some counters were dropped during preprocessing, should use rawSmiles

used after filter_badPatterns.py
'''

def get_smiles_set(file):
    smiles_set=set()
    df = pd.read_csv(file,sep='\t')
    df.dropna(subset=['originalSmiles'],inplace=True)
    for row in df.itertuples():
        mol = Chem.MolFromSmiles(row.originalSmiles)
        if mol:
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
            smiles = Chem.MolToSmiles(mol)
            smiles_set.add(smiles)
    return smiles_set


def get_preprocessed_organismSets(refinedCoconut_df,outputfile,smiles_set):
    preprocessed = refinedCoconut_df[refinedCoconut_df.rawSmiles.isin(smiles_set)]     
    preprocessed.to_csv(outputfile,index=False)
    
    
    
if __name__ == "__main__":
    inputDir = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/origin')
    plants_file = inputDir / 'plants.csv'
    bacteria_file = inputDir / 'bacteria.csv'
    fungi_file = inputDir / 'fungi.csv'
    marine_file = inputDir / 'marine.csv'

    plants_smiles = get_smiles_set(plants_file)
    bacteria_smiles = get_smiles_set(bacteria_file)
    fungi_smiles = get_smiles_set(fungi_file)
    marine_smiles = get_smiles_set(marine_file)

    outputDir = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/noBadPatternSmiles')
    if not outputDir.exists():
        outputDir.mkdir()

    plants_out = outputDir / 'plants_noBadPattern.csv'
    bacteria_out = outputDir / 'bacteria_noBadPattern.csv'
    fungi_out = outputDir / 'fungi_noBadPattern.csv'
    marine_out = outputDir / 'marine_noBadPattern.csv'

    #refinedCoconut
    inputFolder = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/noBadPatternSmiles'
    csv_files = glob.glob(os.path.join(inputFolder, "*.csv"))
    frames = []
    for f in csv_files:
        df_part = pd.read_csv(f,sep='\t')
        frames.append(df_part)
    refinedCoconut_df = pd.concat(frames)
    
    get_preprocessed_organismSets(refinedCoconut_df,plants_out,plants_smiles)
    get_preprocessed_organismSets(refinedCoconut_df,bacteria_out,bacteria_smiles)
    get_preprocessed_organismSets(refinedCoconut_df,fungi_out,fungi_smiles)
    get_preprocessed_organismSets(refinedCoconut_df,marine_out,marine_smiles)
    
    print('Done')
