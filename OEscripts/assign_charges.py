import os
from openeye import oechem
from openeye import oequacpac
import multiprocessing as mp
from datetime import datetime


'''
Script to calculate the partial charges with the AM1BCCELF10 model for the ring systems.
Use after the application of OMEGA.
'''

def assignCharges_LargeFile(inputFile,outputFile):
    jobList = []
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()

    if ifs.open(inputFile):

        for mol in ifs.GetOEMols():
            jobList.append(oechem.OEMol(mol))

    else:
        oechem.OEThrow.Fatal('Unable to open {}.'.format(ifs))

    print('Length of joblist: {}'.format(len(jobList)))

    if ofs.open(outputFile):
        pool = mp.Pool(mp.cpu_count()-8)
        print('Start calculating charges.')
        results = pool.map(assignCharges, jobList)
#         print(results)
        pool.close()
        pool.join()

        for mol in results:
            oechem.OEWriteMolecule(ofs, mol)

    else:
        oechem.OEThrow.Fatal('Unable to create {}.'.format(ofs))

        
def assignCharges_smallFile(inputFile,outputFile):
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()

    if not ifs.open(inputFile):
        oechem.OEThrow.Fatal('Unable to open {}.'.format(ifs))
        
    if not ofs.open(outputFile):
        oechem.OEThrow.Fatal('Unable to create {}.'.format(ifs))
        
    for mol in ifs.GetOEMols():
        mol = assignCharges(oechem.OEMol(mol))
        oechem.OEWriteMolecule(ofs, mol)

        
def assignCharges(mol):
    '''
    Assgins the AM1BCCELF10 partial charges to the given molecule.
    :param mol: An openeye molecule for which the partial charges have to be calculated.
    :return: The openeye molecule with assigned partial charges
    '''
    oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCELF10Charges())
    return mol



if __name__ == "__main__":
    start_time = datetime.now()
    
    inputFolder_npRings = '/data/local/ringsys/202111_3dComparison/assign_charges/queries_input_oeFormat/'
    outputFolder_npRings = '/data/local/ringsys/202111_3dComparison/assign_charges/queries_charges_assigned/'
    if not os.path.exists(outputFolder_npRings):
        os.makedirs(outputFolder_npRings)
    
    jobList = []
    for file in os.listdir(inputFolder_npRings):
        infilePath = inputFolder_npRings+file
        outfilePath = outputFolder_npRings+file.split('.')[0]+'_charges_assigned.oeb.gz'
        jobList.append((infilePath,outfilePath))

    print('Length of joblist: {}'.format(len(jobList)))

    pool = mp.Pool(mp.cpu_count()-8)
    print('Start calculating charges.')
    results = pool.starmap(assignCharges_smallFile, jobList)
    pool.close()
    pool.join()
    
    print('Finished in:')
    print(datetime.now()-start_time)
