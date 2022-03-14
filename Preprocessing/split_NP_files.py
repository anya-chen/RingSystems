import os
from pathlib import Path
import argparse
'''
Helper script to split large files, used before preprocess_molecules.py
'''

def split_large_files(file, filesize):
    '''
    Splits a large file with molecules and ID into smaller files of the given file size and saves them in
    ../smiles_input/databse_n.smi

    :param file: large file
    :param filesize: number of lines in the smaller files
    '''
    outputFolder = Path(file).parent.parent / 'smiles_input'
    if not outputFolder.exists():
        outputFolder.mkdir()
    
    smallfile = None
    with open(file, 'r', encoding='utf8') as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % filesize == 0:
                if smallfile:
                    smallfile.close()
                fileName = os.path.basename(file).split('.')[0]
                outFilename = fileName + "_" + str(lineno//filesize) + ".smi"
                outpath = outputFolder / outFilename
                smallfile = open(outpath, 'w', encoding='utf-8')
                if not lineno == 0:
                    smallfile.write('smiles ID\n')
            smallfile.write(line)
        if smallfile:
            smallfile.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--inputfilepath", type=str, help='Specify input filepath with files to split.',
                   default="/home/ychen/projects/nps_ringsys/20210519_data_prep/NPs/smiles")
    p.add_argument('-n', '--splitnumber', type=int, help='Number of molecules in one file.',
                   default=3000)
    args = p.parse_args()
    for file in os.listdir(args.inputfilepath):
        if file.endswith(".smi"):
            split_large_files(os.path.join(args.inputfilepath, file), args.splitnumber)

