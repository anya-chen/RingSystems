import os
from pathlib import Path
import csv
from preprocess_molecules import *
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
from separate_multicomponents import separate_multicomponents
from filter_badPatterns import *
from remove_duplicate import remove_duplicate_splitted_files


#first convert the sdf to a csv with "smiles id"
drugs = PandasTools.LoadSDF('/data/local/ringsys/202111_approveddrug/approveddrug.sdf')
drugs_csv_path = '/data/local/ringsys/202111_approveddrug/reformated_input/approveddrug.csv'
drugs[['SMILES','DATABASE_ID']].to_csv(drugs_csv_path,index=False,sep=' ')

#preprocessing
preprocess_database(drugs_csv_path)
drugs_dir = '/data/local/ringsys/202111_approveddrug/'
separate_multicomponents(Path(drugs_dir)/'preprocessedSmiles')
badPatterns_filter(Path(drugs_dir)/'separate_multicom')

#remove duplicates
inputFolder = Path(drugs_dir)/'noBadPatternSmiles'
outputFolder = inputFolder.parent / 'uniqueSmiles'
if not outputFolder.exists():
    outputFolder.mkdir()

remove_duplicate_splitted_files(inputFolder,outputFolder,'drug')

