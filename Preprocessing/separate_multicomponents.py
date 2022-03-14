import os
from pathlib import Path
from rdkit import Chem
import filecmp

'''
Script to split multi-component compounds into individual compounds.
Use after preprocess_Zinc.py preprocess_coconut.py
better way: use Chem.GetMolFrags(mol,asMols=True) because there could be '.' but still be one component, 
e.g. CCCC1.C1CCCC, but here there is no such special cases so using '.' is fine.
'''

def separate_multicomponents(inputDir):
    outdir = Path(inputDir).parent / 'separate_multicom'
    if not outdir.exists():
        outdir.mkdir()
    print(inputDir)
    
    for filename in os.listdir(inputDir):
        outfilename = filename.split('.')[0]+'_noMulticom.csv'
        infile = inputDir / filename
        outfile = outdir / outfilename
        with open(infile, 'r') as inputSmilesCSV:
            with open(outfile, 'w') as of:
                header = inputSmilesCSV.readline()
                of.write(header)
                for line in inputSmilesCSV:
                    entry = line.rstrip().split('\t')
                    if len(entry) == 4:
                        preprocessedSmile = entry[1]
                        name = entry[2]
                        smilesList = preprocessedSmile.split('.')
                        smilesList2 = entry[3].split('.')
                        if len(smilesList) ==1:
                            of.write(line)
                        else:
                            for i in range(len(smilesList)):
                                new_name = str(name)+'x'+str(i)
                                of.write(entry[0]+'\t'+smilesList[i]+'\t'+new_name+'\t'+smilesList2[i]+'\n')
            
            

if __name__ == "__main__":
    zinc_processedDir = Path("/home/ychen/projects/nps_ringsys/20210519_data_prep/ZINC/preprocessedSmiles") #no multi-component compounds for zinc
    refinedCoconut_processedDir = Path("/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/preprocessedSmiles")
    removedCoconut_processedDir = Path("/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut/preprocessedSmiles")

    separate_multicomponents(zinc_processedDir)
    separate_multicomponents(refinedCoconut_processedDir)
    separate_multicomponents(removedCoconut_processedDir)
    