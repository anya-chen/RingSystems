import os
from pathlib import Path
from rdkit import Chem
import multiprocessing as mp
from datetime import datetime

'''
This script is used to get smiles that has at least one ring. 
use before get ring systems.
'''

def get_smiles_has_ring(infile,outfile):
    with open(outfile,'w') as of:
        with open(infile,'r') as f:
            header = f.readline()
            of.write(header)
            for line in f:
                entry = line.rstrip().split('\t')
                mol = Chem.MolFromSmiles(entry[0])
                if Chem.GetSSSR(mol) != 0:
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

        outputPath = inputfilePath.parent / 'smilesHasRing'
        if not outputPath.exists():
            outputPath.mkdir()

        for file in os.listdir(inputfilePath):
            filePath = inputfilePath / file
            outfilename = file.split('.')[0] + '_hasRing.csv'
            outfilePath = outputPath / outfilename
            jobList.append((filePath,outfilePath))
    
    pool = mp.Pool(processes=mp.cpu_count() - 2)
    pool.starmap(get_smiles_has_ring, jobList)
    pool.close()
    pool.join()
    
    print('Finished in: ')
    print(datetime.now() - begin_time)
    print('Done')
    