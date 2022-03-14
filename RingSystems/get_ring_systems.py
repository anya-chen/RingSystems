import pandas as pd
from rdkit import Chem
from RingSystemClass import MoleculeWithRings

'''
Script used in get_ZINC_rings.py and get_NP_rings.py
'''

def process_molecules(filepath):
    '''
    for each molecule of the dictionary the ring systems are obtained and stored in a ring system dataframe

    :param moleculeDF: molecule dataframe with ['index','preprocessedSmiles','id','conID','mol_stereoIDs']
    :return: pandas Dataframe for both datasets with the following columns:
    Ring, RingSmiles, RingName, Molecule, MoleculeSmiles, MoleculeName
    '''
    moleculeDF = pd.read_csv(filepath,sep=',')
    moleculeDF = moleculeDF[['preprocessedSmiles','id','conID','mol_stereoIDs']]
    moleculeDF['molecule'] = moleculeDF.apply(lambda row: Chem.MolFromSmiles(row['preprocessedSmiles']),axis=1)
    
    RingsDF = get_ringsDF(moleculeDF)
    RingsDF = remove_duplicates(RingsDF, ['MoleculePreprocessedSmiles', 'RingSmiles'])

    return RingsDF
    

def get_ringsDF(moleculeDF):
    '''
    creates a ring system dataframe of the ring systems of each molecule in moleculeDF

    :param moleculeDF:  dataframe with the columns ''preprocessedSmiles','id','conID','mol_stereoIDs' and 'molecule'
    :return: A ring system dataframe with the columns 'Ring', 'RingSmiles', 'RingName', 'Molecule', 'MoleculeSmiles', 'MoleculeId', 'MoleculeConID' and 'mol_stereoIDs'
    '''
    propertyList = ['Ring', 'RingSmiles', 'RingName', 'Molecule', 'MoleculePreprocessedSmiles', 'MoleculeId','MoleculeConID','mol_stereoIDs']
    ringsDF = pd.DataFrame(columns=propertyList)
    moleculeDF['ringsObject'] = moleculeDF['molecule'].apply(lambda mol: MoleculeWithRings(mol))
    moleculeDF['ringsDict'] = moleculeDF['ringsObject']
    for row in moleculeDF.itertuples():
        ringsDict = row.ringsObject.get_ring_system_dict()
        moleculeName = row.id
        moleculeConID = row.conID
        moleculeStereoID = row.mol_stereoIDs
        moleculeSmiles = row.preprocessedSmiles
        mol = row.molecule
        ringCounter = 0
        for ringMol, ringSmiles in ringsDict.items():
            ringSmiles = ringSmiles
            ringName = str(moleculeName) + "_" + str(ringCounter)
            ringDF = pd.DataFrame([[ringMol, ringSmiles, ringName, mol, moleculeSmiles, moleculeName, moleculeConID,moleculeStereoID]], columns=propertyList)
            ringsDF = ringsDF.append(ringDF, ignore_index=True)

            ringCounter += 1

    return ringsDF


def remove_duplicates(ringsDF, listOfSubset):
    '''
    removes the duplicates by means of the Existing_code columns

    :param ringsDF: dataframe with the columns 'Ring', 'RingSmiles', 'RingName', 'Molecule', 'MoleculeSmiles' and 'MoleculeName'
    :param listOfSubset: columns which are used to remove duplicates
    :return: deduplicated dataframe
    '''
    ringsDF = ringsDF.drop_duplicates(subset=listOfSubset)
    return ringsDF
