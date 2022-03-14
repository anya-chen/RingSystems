import glob
from pathlib import Path
from datetime import datetime
from assign_ids_molecules import *
import multiprocessing as mp
from datetime import datetime



if __name__ == "__main__": 
    begin_time = datetime.now()
    
    folder = '/data/local/ringsys/202111_molecules_noRing'
    
    jobList = []
    for path in os.listdir(folder):
        if path == 'zinc' or path == 'refinedCoconut':
            inputfilePath = Path(folder) / path / 'smilesNoRing'

            outputPath = inputfilePath.parent / 'smilesNoRing_ids'
            if not outputPath.exists():
                outputPath.mkdir()

            for file in os.listdir(inputfilePath):
                filePath = inputfilePath / file
                #if file only have the header should be deleted (or run that before)
                if path == 'refinedCoconut':
                    prefix = str(filePath).split('_')[5]
                else:
                    prefix = str(filePath).split('_')[7]
                outfilename = file.split('.')[0] + '_ids.csv'
                outfilePath = outputPath / outfilename
                if outfilename not in os.listdir(outputPath):
                    jobList.append((filePath,prefix,outfilePath))
    
    #organism Set
    inputFolder = Path('/data/local/ringsys/202111_molecules_noRing/organismSets/smilesNoRing')
    outputFolder = inputFolder.parent / 'smilesNoRing_ids'
    if not outputFolder.exists():
        outputFolder.mkdir()
    
    files = inputFolder.glob('*')
    for file in files:
        outFileName = str(file.stem)+'_ids.csv'
        outputFile = outputFolder / outFileName
        prefix = file.stem.split('_')[3]
        if outfilename not in os.listdir(outputPath):
            jobList.append((file,prefix,outputFile))

    pool = mp.Pool(processes=100)
    pool.starmap(assign_ids, jobList)
    pool.close()
    pool.join()
    print('Finished in: ')
    print(datetime.now() - begin_time)    
    print('Done')