import os
from pathlib import Path
from rdkit import Chem

'''
Script to filter out wrong structures,for all ZINC and Coconut
The input files are tab separated and have the columns rawSmiles, preprocessedSmiles, id, tautomerizedSmiles, with preprocessedSmiles being the SMILES without the tautomer standardization step.
Use after preprocess_zinc.py and preprocess_coconut.py
'''

def badPatterns_filter(inputdir):
    outdir = Path(inputdir).parent / 'noBadPatternSmiles'
    outdir2 = Path(inputdir).parent / 'BadPatternSmiles'
    if not outdir.exists():
        outdir.mkdir()
    if not outdir2.exists():
        outdir2.mkdir()
    print(inputdir)
    
    smarts =['[*;R]=[*;R]=[*;R]','[#6+]']
    patterns = [Chem.MolFromSmarts(sma) for sma in smarts]
    
    badPatternFile = outdir2 / 'badPatternSmiles.csv'
    with open(badPatternFile,'w') as badPatternSmiles:
        for filename in os.listdir(inputdir):
            outfilename = filename.split('.')[0]+'_noBadPattern.csv'
            infile = inputdir / filename
            outfile = outdir / outfilename
            with open(infile, 'r') as inputSmilesCSV:
                with open(outfile, 'w') as of:
                    header = inputSmilesCSV.readline()
                    of.write(header)
                    for line in inputSmilesCSV:
                        entry = line.rstrip().split('\t')
                        if len(entry) == 4:
                            preprocessedSmile = entry[1]
                            mol = Chem.MolFromSmiles(preprocessedSmile)
                            results = map(lambda pattern:mol.HasSubstructMatch(pattern),patterns)
                            if any(list(results)) == False:
                                of.write(line)
                            else:
                                badPatternSmiles.write(line)
    
            print('bad patterns in {} filtered'.format(infile))



if __name__ == "__main__":
    zinc_processedDir = Path("/home/ychen/projects/nps_ringsys/20210519_data_prep/ZINC/separate_multicom")
    refinedCoconut_processedDir = Path("/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/separate_multicom")
    removedCoconut_processedDir = Path("/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut/separate_multicom")

    badPatterns_filter(zinc_processedDir)
    badPatterns_filter(refinedCoconut_processedDir)
    badPatterns_filter(removedCoconut_processedDir)
    
    
    
