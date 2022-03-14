import os
from pathlib import Path
from rdkit import Chem
import multiprocessing as mp
from datetime import datetime

'''
This script is used to get smiles that has no ring.
Use before assign_ids_molecules_noRing.py
'''

def get_smiles_has_noRing(infile,outfile):
    with open(outfile,'w') as of:
        with open(infile,'r') as f:
            header = f.readline()
            of.write(header)
            for line in f:
                entry = line.rstrip().split('\t')
                mol = Chem.MolFromSmiles(entry[0])
                if Chem.GetSSSR(mol) == 0:
                    of.write(line)
    print('finished filtering {}'.format(infile))



if __name__ == "__main__":
    begin_time = datetime.now()
    
    folder = '/home/ychen/projects/nps_ringsys/20210618_preprocessing'
    jobList = []
    for path in os.listdir(folder):
        if path == 'zinc':
            inputfilePath = Path(folder) / path / 'uniqueSmilesCleaned'
#         elif path == 'removedCoconut':
#             pass
        else:
            inputfilePath = Path(folder) / path / 'uniqueSmiles'

        outputPathFolder = Path('/data/local/ringsys/202111_molecules_noRing/') / path 
        if not outputPathFolder.exists():
            outputPathFolder.mkdir()
        outputPath = outputPathFolder / 'smilesNoRing'
        if not outputPath.exists():
            outputPath.mkdir()

        for file in os.listdir(inputfilePath):
            filePath = inputfilePath / file
            outfilename = file.split('.')[0] + '_hasNoRing.csv'
            outfilePath = outputPath / outfilename
            jobList.append((filePath,outfilePath))
    
    pool = mp.Pool(processes=10)
    pool.starmap(get_smiles_has_noRing, jobList)
    pool.close()
    pool.join()
    
    print('Finished in: ')
    print(datetime.now() - begin_time)
    print('Done')
    
