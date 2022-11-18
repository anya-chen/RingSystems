'''
MoleculePreprocessor

.. moduleauthor:: Conrad Stork <stork@zbh.uni-hamburg.de> and Anke Wilm <wilm@zbh.uni-hamburg.de>

'''
import sys
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem import AllChem
from src.neutralizationPatterns import neutralizationReactions
import rdkit.Chem.MolStandardize.rdMolStandardize as rdMolStandardize

class MoleculePreprocessor():
 '''
 MoleculePreprocessor class contains a dictionary where the inital SMILES are mapped to the molecule. The class contains different methods for molecule preparation (e.g. removing salts). All steps are written into a log dict where also problematic compounds are listed.

 :param dictOfMolecules: smiles:mol dict for storage
 :param replaceFilteredSmiles: if smiles are replaced if they are filtered
 :param log: log keys !:invalid molecule, 1:removed, 0:changed, S:salt, W:weight, E:element, N:neutralize, C:canonalize, V:validity check.
 '''
 def __init__(self,dictOfMolecules,log,replaceFilteredSmiles):
  self.dictOfMolecules = dictOfMolecules
  self.replaceFilteredSmiles = replaceFilteredSmiles
  self.log = log

 @classmethod
 def init_with_smiles(cls,listOfSmiles,replaceFilteredSmiles=True,removeSteroInfo=False):
  '''
  Initilizes a MoleculePreprocessor with all SMILES from the given List.

  :param listOfSmiles: list of strings
  :param replaceFilteredSmiles: wheather filterd smiles should be replaced by None
  :return: a initilized MoleculePreprocessor class
  '''
  dictOfMolecules = dict()
  log = {'use':[]}
  for raw_smiles in listOfSmiles:
   if removeSteroInfo:
    try:
     mol = Chem.MolFromSmiles(Chem.CanonSmiles(raw_smiles,useChiral=0))
    except Exception as e:
     if not str(type(e)) == '<class \'Boost.Python.ArgumentError\'>':
      raise e
     mol = None
   else:
    mol = Chem.MolFromSmiles(raw_smiles)
   if not mol:
    log[raw_smiles] = ['!1']
   else:
    log[raw_smiles] = []
   dictOfMolecules[raw_smiles] = mol
  return cls(dictOfMolecules,log,replaceFilteredSmiles)

 def acm_wash(self,minWeight=250,maxWeight=900,allowedAtomNrs={1,5,6,7,8,9,14,15,16,17,34,35,53},tautRemoveStereo=True):
  '''
  Applies the salt_filter, weight_filter, element_filter, canonalize_tautomer and check_molecules_validity

  :param minWeight: lower border for removing cpd
  :param maxWeight: upper border for removing cpd
  :param allowedAtomNrs: Atom number of allowed atoms, cpds containg other than listed ones are removed
  :param tautRemoveStereo: if stereo chemistry should be removed if chiral atoms are inside a tautomeric
         pattern and therefore the stereo isomers are theoretically exchangeable
  '''
  self.salt_filter()
  self.weight_filter(minWeight=minWeight,maxWeight=maxWeight)
  self.element_filter(allowedAtomNrs=allowedAtomNrs)
  self.canonalize_tautomer(removeStereo=tautRemoveStereo)
  self.check_molecules_validity()

 def salt_filter(self):
  '''
  Removes salts according to the rules defined in DOI: 10.1002/cmdc.201700673
  '''
  self.log['use'].append('S')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol: continue
   molFrags = Chem.GetMolFrags(mol,asMols=True)
   smileSplitList = []
   for frag in molFrags:
    smileSplitList.append((frag,frag.GetNumAtoms()))
   sortSmileSplitList = sorted(smileSplitList, key=lambda x:x[1], reverse=True)
   if len(sortSmileSplitList) > 1 and sortSmileSplitList[0][1]*0.7 < sortSmileSplitList[1][1] and (not Chem.MolToSmiles(sortSmileSplitList[0][0],True) == Chem.MolToSmiles(sortSmileSplitList[1][0],True)):
    if self.replaceFilteredSmiles:
     self.dictOfMolecules[key] = None
     self.log[key].append('!1')
    self.log[key].append('S1')
   else:
    if len(sortSmileSplitList) > 1:
     self.log[key].append('S0')
    self.dictOfMolecules[key] = sortSmileSplitList[0][0]
  self.neutralise_charges()

 def weight_filter(self,minWeight=250,maxWeight=900):
  '''
  Sets all molecules with a weight above maxWeight and below minWeight to None.

  :param minWeight: lower border for removing cpd
  :param maxWeight: upper border for removing cpd
  '''
  self.log['use'].append('W')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol:
    continue
   weight = MolWt(mol)
   if weight < minWeight or weight > maxWeight:
    if self.replaceFilteredSmiles:
     self.dictOfMolecules[key] = None
     self.log[key].append('!1')
    self.log[key].append('W1')

 def element_filter(self,allowedAtomNrs={1,5,6,7,8,9,14,15,16,17,34,35,53}):
  '''
  Sets all molecules with elements that are not in allowedAtomNrs to None.

  :param allowedAtomNrs: Atom number of allowed atoms, cpds containg other than listed ones are removed
  '''
  self.log['use'].append('E')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol:
    continue
   atomNrs = set()
   atoms = mol.GetAtoms()
   for atom in atoms:
    atomNrs.add(atom.GetAtomicNum())
   notAllowedAtomsInMolecule = atomNrs - allowedAtomNrs
   if not len(notAllowedAtomsInMolecule) == 0:
    if self.replaceFilteredSmiles:
     self.dictOfMolecules[key] = None
     self.log[key].append('!1')
    self.log[key].append('E1')

 def neutralise_charges(self):
  '''
  Neutralizes the molecules according to the patterns defined in src.ruleSets.neutralizationPatterns.
  '''
  self.log['use'].append('N')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol: continue
   for i,(reactant, product) in enumerate(neutralizationReactions):
    while mol and mol.HasSubstructMatch(reactant):
     rms = AllChem.ReplaceSubstructs(mol,reactant,product)
     mol = rms[0]
     try:
      mol.UpdatePropertyCache()
      if not mol:
       self.log[key].append('N1')
       self.log[key].append('!1')
      self.dictOfMolecules[key] = mol
     except ValueError:
      self.dictOfMolecules[key] = None

 def canonalize_tautomer(self, removeStereo=True):
  '''
  Canonalizes the molecules according to the rules defined in tautomerTransformations. For further details see molvs.
  '''
  taut_params = rdMolStandardize.CleanupParameters()
  if not removeStereo:
    taut_params.tautomerRemoveSp3Stereo = False
    taut_params.tautomerRemoveBondStereo = False
    taut_params.tautomerRemoveIsotopicHs = False
  canon = rdMolStandardize.TautomerEnumerator(taut_params)

  self.log['use'].append('C')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol: continue
   Chem.SanitizeMol(mol)
   molc = canon.Canonicalize(mol)
   if not molc:
    self.log[key].append('C1')
    self.log[key].append('!1')
   self.dictOfMolecules[key] = molc

 def check_molecules_validity(self):
  '''
  Funktion that checks if each rdkit mol is convertable into a smiles and back into a molecule
  '''
  self.log['use'].append('V')
  for key in self.dictOfMolecules:
   molIn = self.dictOfMolecules[key]
   if not molIn: continue
   smi = Chem.MolToSmiles(molIn,True)
   mol = Chem.MolFromSmiles(smi)
   if not mol:
    self.log[key].append('V1')
    self.log[key].append('!1')
    self.dictOfMolecules[key] = mol

 def get_rawsmiles_smiles_dict(self):
  '''
  Returns a dict with rawsmiles mapped to smiles

  :return: dict with rawsmiles mapped to smiles
  '''
  dictOfSmiles = dict()
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol:
    dictOfSmiles[key] = ''
   else:
    dictOfSmiles[key] = Chem.MolToSmiles(mol,True)
  return dictOfSmiles

 def get_rawsmiles_mol_dict(self):
  '''
  Returns a dict with rawsmiles mapped to molecules

  :return: dict with rawsmiles mapped to molecules
  '''
  return self.dictOfMolecules

 def get_log(self):
  '''
  Method that returnes a dict with the raw_smiles as key and a list of error keys

  :return: log dict
  '''
  log = dict(self.log)
  del log['use']
  return log

 def write_log_file(self,filePath=str()):
  '''
  Method that writes the log dict to a file

  :param filePath: path and name of the file to which the log should be written
  '''
  f = open(filePath,'w')
  f.write('log keys !:invalid molecule, 1:removed, 0:changed, S:salt, W:weight, E:element, N:neutralize, C:canonalize, V:validity check.\n')
  f.write('USED FILTERS:'+str(self.log['use']))
  for key in self.log:
   if key == 'use': continue
   f.write(key+' '+str([i for i in self.log[key]])+'\n')
  f.close()

