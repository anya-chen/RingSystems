import os
from openeye import oechem

'''
Script used to concatenate the small molecule and macrocycle files after the application of OMEGA.
'''


def concatenate_files(inputfilepath, outputfile):
    ofs = oechem.oemolostream()
    molList = []
    for file in os.listdir(inputfilepath):
        if file.endswith('oeb.gz'):
            print(file)
            ifs = oechem.oemolistream()
            if ifs.open(os.path.join(inputfilepath, file)):
                for mol in ifs.GetOEMols():
                    molList.append(oechem.OEMol(mol))
            else:
                oechem.OEThrow.Fatal('Unable to open input')

    if ofs.open(outputfile):
        for mol in molList:
            oechem.OEWriteMolecule(ofs, mol)
    else:
        oechem.OEThrow.Fatal('Unable to create output')


inputfilepath = '/home/ychen/projects/nps_ringsys/20210901_3dComparison/assign_charges/zinc_charges_assigned/'
outputfile = '/home/ychen/projects/nps_ringsys/20210901_3dComparison/assign_charges/zinc_charges_assigned/zinc_rings_all.oeb.gz'
concatenate_files(inputfilepath, outputfile)
