from rdkit import Chem

def process_query_ligand(query_arg: str):
    """
    Parses a query ligand from either a SMILES string or a path to an SDF file.
    Canonicalizes the structure to freeze atom indices consistently.

    Note: Hydrogens are kept implicit (not added) as Vina scoring uses heavy atoms only.

    Returns:
        (query_mol, canonical_smiles): RDKit Mol object and canonical SMILES string.
    """
    if query_arg.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(query_arg)
        mol = suppl[0]
        if mol is None:
            raise ValueError(f"Failed to load molecule from {query_arg}")
        smiles = Chem.MolToSmiles(mol)
    else:
        smiles = query_arg

    # Canonicalize systematically to freeze atom indices
    try:
        canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        query_mol = Chem.MolFromSmiles(canonical_smiles)
        # No AddHs() - Vina uses heavy atoms only, H added later if needed
    except Exception as e:
        raise ValueError(f"Failed to parse or canonicalize SMILES '{smiles}': {e}")

    return query_mol, canonical_smiles
