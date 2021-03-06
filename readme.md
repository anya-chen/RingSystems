# RingSystems
## Ring systems in natural products: structural diversity, physicochemical properties, and coverage by synthetic compounds
This is the code used in the cheminformatic analysis for ring systems in natural prodcuts, for more details please see the publication: [Natural Product Reports, 2022, DOI: 10.1039/D2NP00001F](https://doi.org/10.1039/D2NP00001F) 

### Requirements
Anaconda (or minicoda) and Git should be installed.   
Lisence is needed to use [OpenEye](https://www.eyesopen.com) applications and toolkits.  
```
git clone https://github.com/anya-chen/RingSystems  
cd RingSystems  
conda env create -n ringsys -f environment.yml  
conda activate ringsys  
pip install -e .  
```
If you are installing manually/using only certain part:  
- Create ringsys env with python 3.8 and [RDKit](https://www.rdkit.org/)  
```
conda create -n ringsys python=3.8
conda activate ringsys
conda install -c conda-forge rdkit
```  
- Install [ChEMBL Structure Pipline](https://github.com/chembl/ChEMBL_Structure_Pipeline)  
    
```conda install -c conda-forge chembl_structure_pipeline```  

- Install [oepython (openeye toolkits)](https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html)  

```conda install -c openeye openeye-toolkits```  

- Install scikit-learn, numpy, pandas, seaborn...


### Input datasets needed
- <code>Data_prep/get_refined_coconut.py and Data_prep/get_organism_sets.py</code>  
    <code>coconut.sourceNP.csv</code>: from [COCONUT database](https://coconut.naturalproducts.net/download) MongoDB dump, version 2020-10
- <code>Preprocessing/preprocess_Zinc.py</code>  
    <code>zinc20/</code> and <code>zinc_catalogs/</code>: in-stock subset and biogenic sets from [ZINC 20 database](https://zinc20.docking.org/)
- <code>Preprocessing/preprocessing_approveddrug.py</code>  
    <code>approveddrug.sdf</code>: from [DrugBank](https://go.drugbank.com/), version 5.1.8
    
    
### Algorithm to get ring systems from molecules
<code>RingSystems/RingSystemClass.py</code>

Ring systems are defined as all atoms forming a ring, plus any proximate exocyclic atom(s) connected via any type of bond other than a single bond. Two rings sharing at least one atom (i.e. fused and spiro rings) are considered as one ring system.
In order to obtain the ring systems the following algorithm was applied to each chemical structure with at least one ring:
1. Split of molecule into individual rings (with the RDKit function ringInfo). This process results in one or more ring atom sets.
2. If two ring atom sets share at least one atom the sets are fused.
3. The resulting ring systems (i.e. processed ring atom sets) are extended by all atoms directly connected to the ring via any type of bond other than a single bond.
4. All other substituents are replaced by a hydrogen atom.  



### Algorithm to test whether or not two molecules are identical (if there is no evidence that the molecules are not identical)
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
        - Enumerate all possible enantiomers based on all atoms in unspecified_atoms (this results in 2<sup>n</sup> enantiomers, where n is the number of atoms in unspecified_atoms)
        - For each enantiomer:
            - Test whether m1 and m2 can be superposed (with the HasSubstructMatch function in RDKit; this time with useChirality=True) 
                - If yes: 
                    - return TRUE
    - return FALSE


