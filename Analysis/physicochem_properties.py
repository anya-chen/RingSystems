import argparse
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
import pandas as pd
from datetime import datetime

'''
Script to calculate 14 physicochemical properties of the ring systems:
number of nitrogen atoms, number of oxygen atoms, number of chiral centers, molecular weight, number of heavy atoms,
number of hydrogen bond acceptors, number of hydrogen bond donors, logP, topological polar surface area, number of
aromatic atoms, formal charge, number of rings, number of bridgehead atoms, fraction of Csp3 atoms

for np rings:
python physicochem_properties.py
for zinc rings:
python physicochem_properties.py -i /home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/uniqueRingSystems/zinc_uniqueRingSystems_noStereo.smi -np 0
'''

def get_dataframe_from_ring_csv(filepath):
    '''
    Creates a dataframe from the input file

    :param filepath: filepath to input file, a tab separated file with ring system SMILES and ring_conID is expected
    :return: a dataframe
    '''
    ringDF = pd.read_csv(filepath, sep='\t')
    return ringDF


def get_physicochemical_properties(ringDF):
    '''
    Applies all property calculations to the ring systems of the dataframe and stores each property in a new column

    :param ringDF: dataframe with ring systems as SMILES in the column 'ringSmiles_noStereo'
    :return: a dataframe with ring system molecules and their properties
    '''
    PandasTools.AddMoleculeColumnToFrame(ringDF, 'ringSmiles_noStereo', 'ringMolecule')
    ringDF.dropna(inplace=True)
    print('Start calculcating parameters...')
    ringDF['N'] = ringDF['ringMolecule'].apply(get_molecule_composition, args=(7,))
    ringDF['O'] = ringDF['ringMolecule'].apply(get_molecule_composition, args=(8,))
    ringDF['chiral'] = ringDF['ringMolecule'].apply(get_nof_chiral_centers)
    ringDF['MW'] = ringDF['ringMolecule'].apply(get_MW)
    ringDF['heavy_atoms'] = ringDF['ringMolecule'].apply(num_heavy_atoms)
    ringDF['h_acc'] = ringDF['ringMolecule'].apply(num_of_h_acceptors_and_donors, args=(True,))
    ringDF['h_don'] = ringDF['ringMolecule'].apply(num_of_h_acceptors_and_donors, args=(False,))
    ringDF['logP'] = ringDF['ringMolecule'].apply(get_logp)
    ringDF['TPSA'] = ringDF['ringMolecule'].apply(get_TPSA)
    ringDF['numAro'] = ringDF['ringMolecule'].apply(num_aromatic_atoms)
    ringDF['formalCharge'] = ringDF['ringMolecule'].apply(sum_formal_charge)
    ringDF['numRings'] = ringDF['ringMolecule'].apply(num_rings)
    ringDF['bridgeheadAtoms'] = ringDF['ringMolecule'].apply(num_bridgehead_atoms)
    ringDF['frac_csp3'] = ringDF['ringMolecule'].apply(fraction_csp3)

    return ringDF


def get_molecule_composition(mol, requestedAtomicNum):
    '''
    Counts the number of atoms of a given element in the ring system

    :param mol: the ring system molecule
    :param requestedAtomicNum: atomic number of the element for which the occurrence should be counted
    :return: the number of atoms of an element
    '''
    counter = 0
    for atom in mol.GetAtoms():
        atomicNum = atom.GetAtomicNum()
        if atomicNum == requestedAtomicNum:
            counter += 1
    return counter


def get_nof_chiral_centers(mol):
    return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))


def get_MW(mol):
    return round(Descriptors.MolWt(mol), 3)


def num_heavy_atoms(mol):
    return Lipinski.HeavyAtomCount(mol)


def num_of_h_acceptors_and_donors(mol, acc=True):
    if acc:
        return Lipinski.NumHAcceptors(mol)
    else:
        return Lipinski.NumHDonors(mol)


def get_logp(mol):
    return round(Crippen.MolLogP(mol), 3)


def get_TPSA(mol):
    return round(Descriptors.TPSA(mol), 3)


def num_aromatic_atoms(mol):
    numAromaticAtoms = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            numAromaticAtoms += 1
    return numAromaticAtoms


def sum_formal_charge(mol):
    formalCharge = 0
    for atom in mol.GetAtoms():
        formalCharge += atom.GetFormalCharge()
    return formalCharge


def num_rings(mol):
    return Chem.GetSSSR(mol)


def num_bridgehead_atoms(mol):
    return rdMolDescriptors.CalcNumBridgeheadAtoms(mol)


def fraction_csp3(mol):
    return round(Descriptors.FractionCSP3(mol), 3)


def largest_ringsystem_by_column(ringDF, column, outputdir):
    '''
    Writes the SMILES of one ring system with the highest number of a given column in a file with the path
    [outputdir]/[column]_largest_rings.txt

    :param ringDF: dataframe with ring systems and their properties
    :param column: column for which the highest number is requested
    :param outputdir: directory where the result file will be stored
    '''
    maxRow = ringDF[column].argmax()
    print(ringDF.loc[[maxRow]][column])
    print(ringDF.at[maxRow, 'ringSmiles_noStereo'])
    with open(outputdir + column + "_largest_rings.txt", 'a') as of:
        of.write('Database: {}, Ring system: {}\n'.format(ringDF.at[maxRow, 'db'], ringDF.at[maxRow, 'ring_conID']))
        of.write('Smiles max {}: {}\n'.format(column, ringDF.at[maxRow, 'ringSmiles_noStereo']))


def save_rings_df_with_properties(ringDF, outputfilePath):
    '''
    Saves the dataframe with the ring system SMILES and their properties in a CSV file

    :param ringDF: dataframe with ring system SMILES, names, RDKit molecules and their properties
    :param outputfilePath: filepath to result file
    '''
    ringDFWithoutMol = ringDF.drop(columns='ringMolecule')
    ringDFWithoutMol.to_csv(outputfilePath, sep=';', index=False)


def calculate_mean_and_sd(ringDF, outputfilePath):
    '''
    calculates the mean, median, standard deviation and maximum value for each numerical column in the dataframe
    and stores them
    additionally, the portion of ring systems with a MW > 500 DA, without a chiral center, nitrogen and oxygen
    is calculated

    :param ringDF: dataframe with ring system SMILES, names, RDKit molecules and their properties
    :param outputfilePath: filepath to result file
    '''
    with open(outputfilePath, 'w') as of:
        of.write('Property\tMean\tMedian\tStandardDeviation\tMax\n')
        for column in ringDF:
            if ringDF[column].dtype == 'int64' or ringDF[column].dtype == 'float64':
                mean = ringDF[column].mean()
                std = ringDF[column].std()
                median = ringDF[column].median()
                maxidx = ringDF[column].argmax()
                maxVal = ringDF.at[maxidx, column]
                of.write('{}\t{}\t{}\t{}\n'.format(column, mean, median, std, maxVal))
                if column == 'MW':
                    numAbove500 = len(ringDF[ringDF[column] > 500])
                    dfSize = len(ringDF)
                    of.write('{}\t{}\t{}\n'.format(column, '>500', numAbove500/dfSize))
                elif column == 'chiral':
                    numWithoutChiral = len(ringDF[ringDF[column] == 0])
                    dfSize = len(ringDF)
                    of.write('{}\t{}\t{}\n'.format(column, 'w/o_chiral', numWithoutChiral/dfSize))
                elif column == 'N':
                    numNnull = len(ringDF[ringDF[column] == 0])
                    dfSize = len(ringDF)
                    of.write('{}\t{}\t{}\n'.format(column, 'w/o_nitrogen', numNnull / dfSize))
                elif column == 'O':
                    numOnull = len(ringDF[ringDF[column] == 0])
                    dfSize = len(ringDF)
                    of.write('{}\t{}\t{}\n'.format(column, 'w/o_oxygen', numOnull / dfSize))


if __name__ == "__main__":
    start_time = datetime.now()
    
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--inputfile", type=str, help='Specify input file with ring systems',
                   default='/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueRingSystems/np_uniqueRingSystems_noStereo.smi')
    p.add_argument("-o", "--outputdir", type=str, help='Specify output directory for results',
                   default="/home/ychen/projects/nps_ringsys/20210901_analysis/")
    p.add_argument("-np", "--naturalproducts", type=lambda x:bool(int(x)), default=True)

    args = p.parse_args()

    ringsDF = get_dataframe_from_ring_csv(args.inputfile)

    if args.naturalproducts:
        ringsDF.loc[:, 'db'] = 'NP'
    else:
        ringsDF.loc[:, 'db'] = 'ZINC'

    ringsDF = get_physicochemical_properties(ringsDF)

    largest_ringsystem_by_column(ringsDF, 'heavy_atoms', args.outputdir)
    largest_ringsystem_by_column(ringsDF, 'numRings', args.outputdir)

    if args.naturalproducts:
        save_rings_df_with_properties(ringsDF, args.outputdir + 'NP_all_descriptors.csv')
        calculate_mean_and_sd(ringsDF, args.outputdir + 'NP_mean_sd_results.txt')
        print('Finished calculating properties for NP ring systems')
    else:
        save_rings_df_with_properties(ringsDF, args.outputdir + 'ZINC_all_descriptors.csv')
        calculate_mean_and_sd(ringsDF, args.outputdir + 'ZINC_mean_sd_results.txt')
        print('Finished calculating properties for ZINC ring systems')

    print(datetime.now()-start_time)