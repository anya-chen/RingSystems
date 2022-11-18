import chembl_structure_pipeline
from rdkit import RDLogger
from src.MoleculePreprocessor import MoleculePreprocessor

class MoleculePreprocessorExtended(MoleculePreprocessor):
 '''
 Class inheriting from MoleculePreprocessor for the use of the ChEMBL structure pipeline (CSP).
 Merge this class with the MoleculesPreprocessor when web server is updated to newer env. Doing it before would crash the server because env1218 does not have the CSP module installed.
 '''

 def csp_wash(self,minWeight=250,maxWeight=900,allowedAtomNrs={1,5,6,7,8,9,14,15,16,17,34,35,53},
              canonalizeTautomers=True, tautRemoveStereo=True):
  '''
  Applies the weight_filter, element_filter and check_molecules_validity. Further the chembl_structure_pipeline is used to standardize the mols

  :param minWeight: lower border for removing cpd
  :param maxWeight: upper border for removing cpd
  :param allowedAtomNrs: Atom number of allowed atoms, cpds containg other than listed ones are removed
  :param canonalizeTautomers: if tautomers should be transformed in a canonical form
  :param tautRemoveStereo: if stereo chemistry should be removed if chiral atoms are inside a tautomeric
         pattern and therefore the stereo isomers are theoretically exchangeable
  '''
  rdLogger = RDLogger.logger()

  # set c++ log level of rdkit temporarily to ERROR so that csp does not clutter
  # stdout with logging
  rdLogger.setLevel(RDLogger.ERROR)
  self.csp()
  rdLogger.setLevel(RDLogger.INFO)

  self.weight_filter(minWeight=minWeight,maxWeight=maxWeight)
  self.element_filter(allowedAtomNrs=allowedAtomNrs)

  if canonalizeTautomers:
   self.canonalize_tautomer(removeStereo=tautRemoveStereo)

  self.check_molecules_validity()

 def csp(self):
  '''
  Applies the chembl_structure_pipeline (standardizer and get_parent)
  '''
  self.log['use'].append('CSP')
  for key in self.dictOfMolecules:
   mol = self.dictOfMolecules[key]
   if not mol:
    continue
   mol = chembl_structure_pipeline.standardize_mol(mol)
   mol,csvErrorCode = chembl_structure_pipeline.get_parent_mol(mol)
   if csvErrorCode:
    self.log[key].append(csvErrorCode)
   #self.log[key].append('CSP0')
   #TODO: is the replace stuff needed? (line below)
   #if self.replaceFilteredSmiles:
   # self.dictOfMolecules[key] = None
   self.dictOfMolecules[key] = mol


