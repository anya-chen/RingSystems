import glob
from pathlib import Path
from datetime import datetime
from assign_ids_molecules_optimized import *
import multiprocessing as mp


if __name__ == "__main__":
    start_time = datetime.now()
    
    inputFolder = Path('/home/ychen/projects/nps_ringsys/20210618_preprocessing/organismSets/smilesHasRing')
    outputFolder = inputFolder.parent / 'smilesHasRing_ids_optimized'
    if not outputFolder.exists():
        outputFolder.mkdir()
    
    jobList = []
    files = inputFolder.glob('*')
    for file in files:
        print(file)
        outFileName = str(file.stem)+'_ids.csv'
        outputFile = outputFolder / outFileName
        prefix = file.stem.split('_')[3]
        jobList.append((file,prefix,outputFile))
    
    pool = mp.Pool(processes=41)
    pool.starmap(assign_ids, jobList)
    pool.close()
    pool.join()
        
    print(datetime.now()-start_time)
    print('Done')