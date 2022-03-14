import pandas as pd
from rdkit import Chem
from rdkit.Chem import Mol, MolToSmiles, MolFromSmiles, ChiralType, BondStereo, AllChem, CanonSmiles, rdFMCS
import itertools
from pathlib import Path
import os
import multiprocessing as mp
from datetime import datetime


'''
#on going for iteration checking
Script used to assign mol_stereoID for molecules in refinedCoconut and zinc.
    Example: 
    RRSS, mol_stereoID = 1
    RSSS, mol_stereoID = 2
    ??S?, mol_stereoID = [1, 2]

(early stop: 2**27 enumeration of chiral tags and 100 matches, results are the same as no early stop)

use after get_smiles_has_ring.py

mol_conID is already assigned during remove duplicates
'''

def assign_ids(filePath,prefix,outfilePath):
    molDF = pd.read_csv(filePath,sep='\t')
    del molDF['tautomerizedSmiles']
    molDF['molSmiles_noBondStereo'] = molDF.apply(lambda row: row['preprocessedSmiles'].replace('\\','').replace('/',''),axis=1)
    molDF['mol_nFlags'] = molDF.apply(lambda row: count_nFlags_smiles(row.molSmiles_noBondStereo),axis=1)
    groupby_molDF = molDF.groupby('conID')
    
    frames = []
    count = 0
    for name,group in groupby_molDF:
        if len(group) ==1 or set(group.mol_nFlags) == {0}:
            count+=1
            group['mol_stereoIDs'] = prefix+'_'+str(count)
            frames.append(group)
        else:
            group.count().reset_index()
            group.sort_values('mol_nFlags',ascending=False,inplace=True)
            mols = list(group.molSmiles_noBondStereo)
            new_mols = mols
            for i in range(1,len(mols)):
                if len(new_mols)>i:
                    new_mols2 = new_mols[0:i]
                    for j in new_mols[i:]:
                        if not superpose(Chem.MolFromSmiles(new_mols[i-1]),Chem.MolFromSmiles(j)):
                            new_mols2.append(j)
                    new_mols = new_mols2

            count = count + len(new_mols)
            new_group = group[group.molSmiles_noBondStereo.isin(new_mols)]
            new_group.drop_duplicates('molSmiles_noBondStereo',inplace=True)
            counts = list(range(count+1)[-len(new_mols):])
            ids = [prefix+'_'+str(c) for c in counts]
            new_group['mol_stereoIDs'] = ids
            #get dropped back and assign them ids
            dropped = group[~group.id.isin(list(new_group.id))]
            dropped = dropped.apply(lambda row:get_id_for_dropped_mol(row,new_group),axis=1)

            group = pd.concat([new_group,dropped])
            frames.append(group)

    stereoDF_mols = pd.concat(frames)
    stereoDF_mols.to_csv(outfilePath,columns=['preprocessedSmiles','id','conID','mol_stereoIDs'])
    
def count_nFlags_smiles(smi):
    #annotate how many chiral tags in the smiles
    smi = smi.replace('@@','@')
    nFlags = smi.count('@')
    return nFlags


def superpose(m1, m2):
    '''
    this algorithm works well if both molecules m1 and m2 have a large number
    of unspecified stereo centers
    not consider bond stereochemistry, so the input molecules should have all removed bond stereo
    '''
  # quick exit: are both molecules the same if we remove all stereo information?
    if MolToSmiles(m1, isomericSmiles=False) != MolToSmiles(m2, isomericSmiles=False):
        return False
    
    #another quick exit: any of the two molecules has 0 chiral tags
    if count_nFlags_smiles(MolToSmiles(m1)) == 0 or count_nFlags_smiles(MolToSmiles(m2)) == 0:
        return True
    
    b = 0
    for match in m1.GetSubstructMatches(m2, uniquify=False, useChirality=False):
        b +=1
        m1_copy = Mol(m1)
        m2_copy = Mol(m2)
        
        # gather all atom pairs where only one atom has a specified stereo center
        relevant_atom_pairs = [(a1, a2) for a1, a2 in zip([m1_copy.GetAtomWithIdx(idx) for idx in match],
                                                          m2_copy.GetAtoms())
                               if (a1.GetChiralTag() == ChiralType.CHI_UNSPECIFIED or
                                   a2.GetChiralTag() == ChiralType.CHI_UNSPECIFIED) and
                               a1.GetChiralTag() != a2.GetChiralTag()]

        # make sure that each tuple starts with the unspecified atom
        relevant_atom_pairs = [(a1, a2) if a1.GetChiralTag() == ChiralType.CHI_UNSPECIFIED
                               else (a2, a1)
                               for (a1, a2) in relevant_atom_pairs]

        # we want to enforce that m1 and m2 have the same CIP label on the stereo centers
        # unfortunately, we don't know which chiral tag (@, @@) corresponds to which CIP label
        # for that reason, we have to try all combinations
        a = 0
        for chiral_tag_combination in itertools.product((ChiralType.CHI_TETRAHEDRAL_CW,
                                                         ChiralType.CHI_TETRAHEDRAL_CCW),
                                                        repeat=len(relevant_atom_pairs)):

            # assign chiral tag combination to atoms
            for (a1, a2), flag in zip(relevant_atom_pairs, chiral_tag_combination):
                # we made sure that the first atom is unspecified
                a1.SetChiralTag(flag)

            if m1_copy.HasSubstructMatch(m2_copy, useChirality=True):
                return True

            if a == 2**27:  #2**27, ~6.5h, early stop after 2**27 enumerations, 2**20 won't capture a few
                return False
            a += 1
            
        if b == 100: #early stop after checking 100 matches; 1000 needs almost 30mins
            return False
        
    return False
 
    
def get_id_for_dropped_mol(row,new_group):
    row['mol_stereoIDs'] = []
    for r in new_group.itertuples():
        if superpose(Chem.MolFromSmiles(r.molSmiles_noBondStereo),Chem.MolFromSmiles(row.molSmiles_noBondStereo)):
            row['mol_stereoIDs'].append(r.mol_stereoIDs)
    return row


    
if __name__ == "__main__": 
    begin_time = datetime.now()
    
    folder = '/home/ychen/projects/nps_ringsys/20210618_preprocessing'
    
    jobList = []
    for path in os.listdir(folder):
        if path == 'zinc' or path == 'refinedCoconut':
            inputfilePath = Path(folder) / path / 'smilesHasRing'

            outputPath = inputfilePath.parent / 'smilesHasRing_ids_optimized'
            if not outputPath.exists():
                outputPath.mkdir()

            for file in os.listdir(inputfilePath):
                filePath = inputfilePath / file
                if path == 'refinedCoconut':
                    prefix = str(filePath).split('_')[5]
                else:
                    prefix = str(filePath).split('_')[7]
                outfilename = file.split('.')[0] + '_ids.csv'
                outfilePath = outputPath / outfilename
                if outfilename not in os.listdir(outputPath):
                    jobList.append((filePath,prefix,outfilePath))

    pool = mp.Pool(processes=mp.cpu_count() - 2)
    pool.starmap(assign_ids, jobList)
    pool.close()
    pool.join()
    print(datetime.now() - begin_time)
    print('Done')
