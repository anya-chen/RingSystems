import os
import glob
from pathlib import Path
import pandas as pd
import csv
from preprocess_molecules import preprocess_database
import multiprocessing as mp

'''
The script to process ZINC, including remove IDs in biogenic subset, common structure preprocessing and tautomerization in order to remove potential natural products from the ZINC dataset. 
Use before remove_NPs_from_ZINC.py
'''


if __name__ == "__main__":

    inputFolder = Path("/home/ychen/projects/nps_ringsys/20210420_original_databases/zinc20")
    outputFolder = Path('/home/ychen/projects/nps_ringsys/20210519_data_prep/ZINC')
    if not outputFolder.exists():
        outputFolder.mkdir()
        
    #Remove NP related IDs (by biogenic subset IDs and biogenic suppliers IDs from ZINC)
    biogenic_inputFolder = Path("/home/ychen/projects/nps_ringsys/20210420_original_databases/zinc_catalogs/")
    bioge = pd.read_table(biogenic_inputFolder / 'biogenic/biogenic.smi', sep=' ', names=['smiles','id'])
    bioge_idSet = set(bioge.id.tolist())
    
    bioge_supplier_idSet = set()
    suplierFolder = biogenic_inputFolder / 'biogenic_suppliers'
    for file in suplierFolder.iterdir():
        dataset = pd.read_table(file,sep=' ')
        ids = set(dataset.zinc_id.tolist())
        bioge_supplier_idSet |= ids 
    
    np_related_idSet = bioge_idSet | bioge_supplier_idSet
    print("There are {} ZINC IDs from biogenic subset and biogenic suppliers' set".format(len(np_related_idSet)))
    
    
    filteredFolder = outputFolder / 'smiles_input'
    if not filteredFolder.exists():
        filteredFolder.mkdir()
    
    for file in os.listdir(inputFolder):
        if file.endswith('.smi'):
            zincFile = inputFolder/file
            outFile = file.split('.')[0]+'_filteredIDs.smi'
            filteredFile = filteredFolder / outFile
            with open(filteredFile,'w',encoding='utf8') as f:
                with open(zincFile, 'r', encoding='utf8') as smilesFile:
                    smilesCsv = csv.reader(smilesFile,delimiter=' ')
                    for row in smilesCsv:
                        if not row[1] in np_related_idSet:
                            f.write('{} {}\n'.format(row[0],row[1]))
                smilesFile.close()
            f.close()
            print('successfully remove NP related IDs and save filtered file {}'.format(filteredFile))

            
    #structure preprocess
    preprocessedDir = outputFolder / 'preprocessedSmiles'
    logdir = outputFolder / 'log'
    if not preprocessedDir.exists():
        preprocessedDir.mkdir()
    if not logdir.exists():
        logdir.mkdir()
    
    databaseFilesJobList = []
    for file in os.listdir(filteredFolder):
        databaseFilesJobList.append(os.path.join(filteredFolder, file))
    pool = mp.Pool(processes=mp.cpu_count()-1)
    pool.map(preprocess_database, databaseFilesJobList)
    pool.close()
    pool.join()
    print('Done')
