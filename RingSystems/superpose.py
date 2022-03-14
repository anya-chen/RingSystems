from rdkit import Chem
from rdkit.Chem import Mol, MolToSmiles, MolFromSmiles, ChiralType, BondStereo, AllChem, CanonSmiles, rdFMCS
import itertools

def superpose(m1, m2):
    '''
    this algorithm works well if both molecules m1 and m2 have a large number
    of unspecified stereo centers
    don't consider bond stereochemistry, so the input molecules should have all removed bond stereo
    '''
  # quick exit: are both molecules the same if we remove all stereo information?
    if MolToSmiles(m1, isomericSmiles=False) != MolToSmiles(m2, isomericSmiles=False):
        return False

    #another quick exit: any of the two molecules has 0 chiral tags
    if count_nFlags_smiles(MolToSmiles(m1)) == 0 or count_nFlags_smiles(MolToSmiles(m2)) == 0:
        return True
    
    for match in m1.GetSubstructMatches(m2, uniquify=False, useChirality=False):
        m1_copy = Mol(m1)
        m2_copy = Mol(m2)
        
        # gather all atom pairs where only one atom has a specified stereo center
        relevant_atom_pairs = [(a1, a2) for a1, a2 in zip([m1_copy.GetAtomWithIdx(idx) for idx in match],
                                                          m2_copy.GetAtoms())
                               if (a1.GetChiralTag() == ChiralType.CHI_UNSPECIFIED or
                                   a2.GetChiralTag() == ChiralType.CHI_UNSPECIFIED) and
                               a1.GetChiralTag() != a2.GetChiralTag()]

        # make sure that each tuple starts with the unspecified atom
        relevant_atom_pairs = [(a1, a2) if a1.GetChiralTag() == ChiralType.CHI_UNSPECIFIED
                               else (a2, a1)
                               for (a1, a2) in relevant_atom_pairs]

        # we want to enforce that m1 and m2 have the same CIP label on the stereo centers
        # unfortunately, we don't know which chiral tag (@, @@) corresponds to which CIP label
        # for that reason, we have to try all combinations
        for chiral_tag_combination in itertools.product((ChiralType.CHI_TETRAHEDRAL_CW,
                                                         ChiralType.CHI_TETRAHEDRAL_CCW),
                                                        repeat=len(relevant_atom_pairs)):

                # assign chiral tag combination to atoms
                for (a1, a2), flag in zip(relevant_atom_pairs, chiral_tag_combination):
                    # we made sure that the first atom is unspecified
                    a1.SetChiralTag(flag)

                if m1_copy.HasSubstructMatch(m2_copy, useChirality=True):
                    return True
    return False

