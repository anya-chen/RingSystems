from rdkit import Chem


def save_ring_systems(ringsDF, filename):
    '''
    Used in get_NP_rings.py and get_ZINC_rings.py. Saves the ring systems in three different formats:
    - only the ring systems as SMILES and their names in a .smi file
    - the original molecule SMILES and names and the ring systems SMILES and names in a .csv file
    - the ring systems in a .sdf.

    :param ringsDF: dataframe with ring systems and the columns 'Ring', 'RingSmiles', 'RingName', 'Molecule',
    'MoleculeSmiles' and 'MoleculeName'
    :param filename: the result file name without the file extension
    '''
    save_rings_as_smiles(ringsDF, filename + '.smi')
    save_rings_with_molecules_as_smiles(ringsDF, filename + '.csv')
    save_rings_as_sdf(ringsDF, filename + ".sdf")
    

def save_rings_as_smiles(ringsDF, filePath):
    '''
    Saves the SMILES and the ID of the ring systems in a space separated file with .smi ending

    :param ringsDF: dataframe with the rdkit molecule, name and SMILES of the origin molecule and the ring systems
    :param filePath: path to file where the ring system should be stored
    '''
    ringsDF.to_csv(filePath, sep=' ', index=False, columns=['RingSmiles','RingName'])


def save_rings_with_molecules_as_smiles(ringsDF, filePath):
    '''
    Saves the ring systems with their origin molecule in a tab separated file with .csv ending in the order
    ringSmiles, ringName, moleculeSmiles, moleculeName

    :param ringsDF: dataframe with the rdkit molecule, name and SMILES of the origin molecule and the ring systems
    :param filePath: path to file where the ring system should be stored
    '''
    ringsDFwithoutMols = ringsDF.drop(columns=['Molecule', 'Ring'])
    ringsDFwithoutMols.to_csv(filePath, sep='\t', index=False)


def save_rings_as_sdf(ringsDF, filename):
    '''
    Saves the ring systems in a 2D .sdf

    :param ringsDF: dataframe with the rdkit molecule, name and SMILES of the origin molecule and the ring systems
    :param filename: path to file where the ring system should be stored
    '''
    failcounter = 0
    failedMoleculesSDF = Chem.SDWriter('failedMolecules.sdf')
    failedMoleculesSmilesFile = 'failedMolecules.smi'
    sdFile = Chem.SDWriter(filename)

    for idx in ringsDF.index:
        ring = ringsDF.at[idx, 'Ring']
        ringName = ringsDF.at[idx, 'RingName']
        try:
            ring.SetProp("_Name", ringName)
            sdFile.write(ring)
        except Exception as e:
            print(e)
            failcounter += 1
            failedMoleculesSDF.write(ringsDF['Molecule'][idx])
            smilesFile = open(failedMoleculesSmilesFile, 'a')
            smilesFile.write("{},{}\n".format(ringsDF.at[idx, 'MoleculePreprocessedSmiles'], ringsDF.at[idx, 'RingSmiles']))
            smilesFile.close()

    print("fails =", failcounter)

