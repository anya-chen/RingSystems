import os
from pathlib import Path
from assign_ids_molecules import assign_ids
from get_NP_rings import get_NP_rings
from get_ring_conID_split import assign_conID_split_files
from assign_ring_stereoIDs import assign_ring_stereoIDs
from get_unique_rings_strict_count import *


drugs_dir = '/data/local/ringsys/202111_approveddrug/'


#assgin id for molecules
infilePath = Path(drugs_dir) / 'uniqueSmiles'
outputPath = infilePath.parent / 'smiles_ids'
if not outputPath.exists():
    outputPath.mkdir()

for file in os.listdir(infilePath):
    filePath = infilePath / file
    print(filePath)
    prefix = str(filePath).split('_')[2].split('.')[0]
    outfilename = file.split('.')[0] + '_ids.csv'
    outfilePath = outputPath / outfilename
    assign_ids(filePath,prefix,outfilePath)


#get ring systems
inputPath = Path(drugs_dir) / 'smiles_ids'
saveRingsPath = Path(inputPath).parent / 'ringSystems'
if not saveRingsPath.exists():
    saveRingsPath.mkdir()
get_NP_rings(inputPath, saveRingsPath)


# assign ring conID
inputFolder = Path(drugs_dir) /'ringSystems'
outputFolder = inputFolder.parent / 'ringSystems_conID'
if not outputFolder.exists():
    outputFolder.mkdir()

assign_conID_split_files(inputFolder,outputFolder,'drug',size=300)


# assign ring stereoIDs
inputfilePath = Path(drugs_dir) / 'ringSystems_conID'

outputPath = inputfilePath.parent / 'ringSystems_stereoID'
if not outputPath.exists():
    outputPath.mkdir()

for file in os.listdir(inputfilePath):
    filePath = inputfilePath / file
    print(filePath)
    prefix = str(filePath).split('_')[4].split('.')[0]
    outfilename = file.split('.')[0] + '_ids.csv'
    outfilePath = outputPath / outfilename
    assign_ring_stereoIDs(filePath,prefix,outfilePath)


# get unique rings
drugIdrugutDir = Path(drugs_dir) / 'ringSystems_stereoID'
drugOutputDir  = Path(drugs_dir) / 'uniqueRingSystems'
if not drugOutputDir.exists():
    drugOutputDir.mkdir()
drugOutput1 = drugOutputDir / 'drug_uniqueRingSystems_noStereo.csv'
drugOutput2 = drugOutputDir / 'drug_uniqueRingSystems_Stereo.csv'

get_unique_ringsys_noStereo(drugIdrugutDir,drugOutput1,'drug')
get_unique_ringsys_Stereo(drugIdrugutDir,drugOutput2,'drug')

