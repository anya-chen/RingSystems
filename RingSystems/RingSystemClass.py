from rdkit import Chem
from rdkit.Chem import AllChem


class MoleculeWithRings:
    '''
    MoleculeWithRings class contains a molecule and the corresponding ring systems
    workflow:
    extract ringatom sets with atom ids from molecule with ringInfo (_get_list_of_ring_atom_sets)
    extend ringatom sets with atoms which are directly connected to the ring via other than a single bond
    (_extend_ringSet)
    search for bonds that have one atom in a ring and one outlying -> breakpoints for fragments
    fragment molecule on breakpoints (_get_fragments_from_mol)
    Replace dummy atoms of the fragmentation step by hydrogen atoms (_get_ring_mol_without_dummy_atom)

    Used in get_ring_systems.py
    '''
    def __init__(self, mol):
        '''
        :param mol: molecule for which the ring systems should be created.
        '''
        self.mol = mol
        self.ringList = self._get_ring_system_list()

    def get_ring_system_dict(self):
        '''
        creates a dictionary which contains the ring system molecules and their SMILES

        :return: ring system molecule:SMILES dictionary
        '''
        ringDict = {}
        for ringSystem in self.ringList:
            ringDict[ringSystem] = self._get_smiles_of_ring_system_without_dummy_atoms(ringSystem)  
        return ringDict

    def _get_ring_system_list(self):
        '''
        Method that takes the rings of a molecule and returns a list of the ring systems

        :return: list of ring system molecules
        '''
        ringAtomsSets = self._get_list_of_ring_atom_sets()
        for ringAtomsSet in ringAtomsSets:
            ringAtomsSet = self._extend_ringSet(ringAtomsSet)
        ringMolList = self._get_fragments_from_mol(ringAtomsSets)
        ringMolListWithoutDummyAtom = self._get_ring_mol_without_dummy_atom(ringMolList)
        return ringMolListWithoutDummyAtom

    def _get_list_of_ring_atom_sets(self):
        '''
        Recieves rings with mol.GetRingInfo() and extends the single ring atom sets to ring system atoms sets if
        two rings share at least one atom

        :return: a list of ring system atom sets
        '''
        ringInfo = self.mol.GetRingInfo()
        listOfRingAtomSets = []
        for atomRing in ringInfo.AtomRings():
            ringAtomSet = set(atomRing)
            newList = []
            for atomSet in listOfRingAtomSets:
                BridgeAtomSet = ringAtomSet.intersection(atomSet)
                if len(BridgeAtomSet):
                    ringAtomSet = (ringAtomSet.union(atomSet))
                else:
                    newList.append(atomSet)
            newList.append(ringAtomSet)
            listOfRingAtomSets = newList
        return listOfRingAtomSets

    def _extend_ringSet(self, ringAtomsSet):
        '''
        include all proximate atoms connected via any bond other than a single bond for ringAtomsSet

        :return: extended ringAtomsSet
        '''
        extendSet = set()
        for moleculeAtom in self.mol.GetAtoms():
            atomIdx = moleculeAtom.GetIdx()
            if atomIdx in ringAtomsSet:
                continue
            else:
                neighborAtoms = moleculeAtom.GetNeighbors()
            for neighbor in neighborAtoms:
                if (neighbor.GetIdx() in ringAtomsSet) and neighbor.IsInRing():
                    if self._is_not_single_bond(moleculeAtom, neighbor):
                        extendSet.add(moleculeAtom.GetIdx())
        ringAtomsSet.update(extendSet)
        return ringAtomsSet

    def _is_not_single_bond(self, atom, neighbor):
        atomIdx = atom.GetIdx()
        neighborIdx = neighbor.GetIdx()
        bond = self.mol.GetBondBetweenAtoms(atomIdx, neighborIdx)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return True
        else:
            return False

    def _get_fragments_from_mol(self, fragmentAtomSetsList):
        '''
        get fragments indicated by atom index sets as rdkit mol objects using the method fragmentOnBonds
        note: fragment by atom index set results in missing bonds for some atoms (e.g. aromatic N) therefore I'm using fragmentOnBonds

        :param fragmentAtomSetsList: list with fragments as atom index sets
        :return: Rdkit molecules of ring systems with any atoms at breakpoints
        '''
        ringMolList = []
        for fragmentSet in fragmentAtomSetsList:
            breakpoints = set()
            markingAtomInRing = self.mol.GetAtomWithIdx(list(fragmentSet)[0])  # choose atom in ring system to mark the current ring system
            Chem.SetSupplementalSmilesLabel(markingAtomInRing, '<xxx>')
            for bond in self.mol.GetBonds():
                begin = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                breakpointIdx = self._is_bond_breakpoint(bond, begin, end,
                                                   fragmentSet)
                if breakpointIdx is not None:
                    breakpoints.add(breakpointIdx)
            if breakpoints:  # molecule has to be fragmented
                breakpoints = list(breakpoints)
                fragments = Chem.FragmentOnBonds(self.mol, breakpoints)
                Chem.SanitizeMol(fragments, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                ringFragmentList = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=False)  # sanitizeFrags has to be False because of incorrect aromaticity
                for molFragment in ringFragmentList:
                    for atom in molFragment.GetAtoms():
                        if atom.GetIsAromatic() and not atom.IsInRing():
                            atom.SetIsAromatic(False)
                        if Chem.GetSupplementalSmilesLabel(atom) == '<xxx>':
                            Chem.SetSupplementalSmilesLabel(atom, "")
                            Chem.SetSupplementalSmilesLabel(markingAtomInRing, '')  # remove tag
                            ringMolList.append(molFragment)
            else:  # entire molecule consists of ring system
                Chem.SetSupplementalSmilesLabel(markingAtomInRing, '')  # remove tag
                ringMolList.append(self.mol)
                return ringMolList
        return ringMolList

    @staticmethod
    def _is_bond_breakpoint(bond, beginAtom, endAtom, fragmentSet):
        if (beginAtom in fragmentSet and endAtom not in fragmentSet) or (
                beginAtom not in fragmentSet and endAtom in fragmentSet):
            return bond.GetIdx()

    def _get_ring_mol_without_dummy_atom(self, ringMolList): 
        '''
        replace any atoms by hydrogen and remove aromatic tag if necessary

        :return: ring system without any atoms
        '''
        ringMolListWithoutDummyAtom = list()
        for ring in ringMolList:
            ringWithoutDummyAtom = AllChem.ReplaceSubstructs(ring,
                                                         Chem.MolFromSmiles('*'),
                                                         Chem.MolFromSmiles('[H]'),
                                                         True)
            ringWithoutDummyAtom = ringWithoutDummyAtom[0]

            # if two rings are connected by a double bond only the immediate atom is considered and can be marked as aromatic
            # which results in an invalid molecule. Therefore the aromatic tag has to be removed.

            for atom in ringWithoutDummyAtom.GetAtoms():
                if atom.GetIsAromatic() and not atom.IsInRing():
                    atom.SetIsAromatic(False)
                    neighbors = atom.GetNeighbors()
                    for neighbor in neighbors:
                        bond = ringWithoutDummyAtom.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetIsAromatic():
                            bond.SetIsAromatic(False) 
            smiles = Chem.MolToSmiles(ringWithoutDummyAtom)  #To keep the ring systems the same as in the molecule, don't convert mol-smi-mol.
            ringWithoutDummyAtom = Chem.MolFromSmiles(smiles) 

            ringMolListWithoutDummyAtom.append(ringWithoutDummyAtom)
        return ringMolListWithoutDummyAtom
    

    @staticmethod
    def _get_smiles_of_ring_system_without_dummy_atoms(ringSystem):
        '''
        Method returns SMILES after the replacement of dummy atoms

        :return: SMILES of ring system
        '''
        if ringSystem is not None:
            smilesWithoutDummyAtoms = Chem.MolToSmiles(ringSystem)
            return smilesWithoutDummyAtoms
        return None