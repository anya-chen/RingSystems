# RingSystems
## Ring systems in natural products: structural diversity, physicochemical properties, and coverage by synthetic compounds
This is the code used in the cheminformatics analysis for ring systems in natural prodcuts, for more details please see the publication: XXX (to be added) 

### Requirements
Anaconda (or minicoda) and Git should be installed. 
Lisence is needed to use [OpenEye](https://www.eyesopen.com) applications and toolkits.
<code>git clone https://github.com/anya-chen/RingSystems
cd RingSystems
conda env create -n ringsys -f environment.yml
conda activate ringsys
pip install -e .
</code>
If you are installing manually:
Create ringsys env with python 3.8 and rdkit
<code>conda create -n ringsys python=3.8
conda activate ringsys
conda install -c conda-forge rdkit</code>

Install chembl_structure_pipline
<code>conda install -c conda-forge chembl_structure_pipeline</code>
Install oepython
<code>conda install -c openeye openeye-toolkits</code>
Install scikit-learn, numpy, pandas, seaborn


### Input datasets needed
- Data_prep/get_refined_coconut.py and Data_prep/get_organism_sets.py  
    coconut.sourceNP.csv: from https://coconut.naturalproducts.net/download MongoDB dump, version 2020-10
- Preprocessing/preprocess_Zinc.py  
    zinc20/ and zinc_catalogs/: in-stock subset and biogenic sets from ZINC 20 database: https://zinc20.docking.org/
- Preprocessing/preprocessing_approveddrug.py  
    approveddrug.sdf: from https://go.drugbank.com/, version 5.1.8
    
    
### Algorithm to get ring systems from molecules
<code>RingSystems/RingSystemClass.py</code>

Ring systems are defined as all atoms forming a ring, plus any proximate exocyclic atom(s) connected via any type of bond other than a single bond. Two rings sharing at least one atom (i.e. fused and spiro rings) are considered as one ring system.
In order to obtain the ring systems the following algorithm was applied to each chemical structure with at least one ring:
1. Split of molecule into individual rings (with the RDKit function ringInfo). This process results in one or more ring atom sets.
2. If two ring atom sets share at least one atom the sets are fused.
3. The resulting ring systems (i.e. processed ring atom sets) are extended by all atoms directly connected to the ring via any type of bond other than a single bond.
4. All other substituents are replaced by a hydrogen atom.
From a single compound more than one ring system may be derived. Multiple occurrences of a specific ring system in one and the same molecule increase the count of ring systems by 1 



### Algorithm to identify if two molecules are identical (if there is no evidence that the molecules are not identical)
<code>RingSystems/superpose.py</code>

In the scenario considering stereochemistry (i.e. tetrahedral atom configuration), pairs of molecules were tested for identity according to a procedure that builds on the evidence-based approach. The procedure returns TRUE for a pair of molecules, m1 and m2, if the two molecules are identical (more accurately, if there is no evidence that the molecules are not identical):
- If the constitution of m1 and m2 is distinct (based on their SMILES notations, with any stereochemical information removed):
    - return FALSE
- If the constitution of m1 and m2 is identical (based on their SMILES notations, with any stereochemical information removed):
    - Generate all possible substructure matches between m1 and m2 (using the GetSubstructMatches function of RDKit; stereochemical information disregarded with useChirality=False)
    - For each substructure match:
        - For each pair of matching atoms:
            - If the configuration of exactly one atom is not specified:
                - Add the unspecified atom to unspecified_atoms (a list of atoms for which their configuration will be enumerated)
        - Enumerate all possible enantiomers based on all atoms in unspecified_atoms (this results in 2n enantiomers, where n is the number of atoms in unspecified_atoms)
        - For each enantiomer:
            - Test whether m1 and m2 can be superposed (with the HasSubstructMatch function in RDKit; this time with useChirality=True) 
                - If yes: 
                    - return TRUE
    - return FALSE


