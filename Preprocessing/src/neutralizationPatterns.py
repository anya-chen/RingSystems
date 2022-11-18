from rdkit import Chem

patts= (
    # Imidazoles
    ('[n+;H]','n'),
    # Amines
    ('[N+;!H0]','N'),
    # Carboxylic acids and alcohols
    ('[$([O-]);!$([O-][#7])]','O'),
    # Thiols
    ('[S-;X1]','S'),
    # Sulfonamides
    ('[$([N-;X2]S(=O)=O)]','N'),
    # Enamines
    ('[$([N-;X2][C,N]=C)]','N'),
    # Tetrazoles
    ('[n-]','[nH]'),
    # Sulfoxides
    ('[$([S-]=O)]','S'),
    # Amides
    ('[$([N-]C=O)]','N'),
    )
neutralizationReactions = [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]
