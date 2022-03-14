import os
from pathlib import Path
from datetime import datetime

'''
split 3 large files into smaller parts of one file has 50k rows
'''

def split_large_file(infile,filesize =50000):
    outputFolder = Path(infile).parent
    smallfile = None
    with open(infile,'r') as bigfile:
        header = bigfile.readline()
        print(header)
        for lineno, line in enumerate(bigfile):
            if lineno % filesize == 0:
                if smallfile:
                    smallfile.close()
                fileName = os.path.basename(infile).split('.')[0]
                outFilename = fileName + '_' + str(lineno//filesize) +'.csv'
                outpath = outputFolder /outFilename
                smallfile = open(outpath,'w')
                smallfile.write(header)
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    
    bigfileFolder = Path(infile).parent.parent / 'uniqueSmiles_largeFiles'
    if not bigfileFolder.exists():
        bigfileFolder.mkdir()
    os.system('mv ' + str(infile) +' '+ str(bigfileFolder))


if __name__ == "__main__":
    begin_time = datetime.now()
    refinedCoconutLarge = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/refinedCoconut/uniqueSmiles/refinedCoconut_noVariance_uniq.csv'
    split_large_file(refinedCoconutLarge)

    removedCoconutLarge = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/removedCoconut/uniqueSmiles/removedCoconut_noVariance_uniq.csv'
    split_large_file(removedCoconutLarge)

    zincLarge = '/home/ychen/projects/nps_ringsys/20210618_preprocessing/zinc/uniqueSmilesCleaned/NP_removed_zinc_noVariance_uniq.csv'
    split_large_file(zincLarge)
    
    print('Finished in: ')
    print(datetime.now() - begin_time)