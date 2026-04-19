process SplitLigands {

    tag "${ligands}:batch=${batch_size}"

    publishDir(
        "${params.outputs}/ligands/splits", 
        mode: 'copy', 
        saveAs: { v -> "${ligands.simpleName}-${v}"},
    )

    input:
    path ligands
    val batch_size

    output:
    path "lib_*.{sdf,smi}"

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    import gzip
    from itertools import batched

    from rdkit import Chem

    if "${ligands}".endswith(".gz"):
        opener = partial(gzip.open, mode="rb")
    else:
        opener = partial(open, mode="rb")

    if "${ligands}".endswith((".sdf", ".sdf.gz")):
        supplier = partial(Chem.ForwardSDMolSupplier, removeHs=False)
    elif "${ligands}".endswith((".smi", ".smi.gz")):
        supplier = partial(Chem.SmilesMolSupplierFromText, titleLine=False)
    else:
        raise ValueError("Filename '${ligands}' should end with .smi or .sdf (± .gz)")

    with opener("${ligands}") as f_in:
        if "${ligands}".endswith((".sdf", ".sdf.gz")):
            supp_input = f_in
        elif "${ligands}".endswith((".smi", ".smi.gz")):
            supp_input = "\\n".join(line.decode("utf-8") for line in f_in)
        else:
            raise ValueError("Filename '${ligands}' should end with .smi or .sdf (± .gz)")
        with supplier(supp_input) as supp:
            for i, mols in enumerate(batched(supp, ${batch_size})):
                with open(f"lib_{i:06d}.sdf", "a") as f:
                    for j, mol in enumerate(mols):
                        if mol is None:
                            continue
                        if mol.HasProp("zinc_id"):
                            this_name = mol.GetProp("zinc_id")
                        elif mol.HasProp("name"):
                            this_name = mol.GetProp("name")
                        elif mol.HasProp("_Name"):
                            this_name = mol.GetProp("_Name")
                        else:
                            this_name = f"${ligands.simpleName}-{(i * ${batch_size} + j):08d}"
                        mol.SetProp("_Name", this_name)
                        mol.SetProp("name", this_name)
                        if not mol.HasProp("smiles"):
                            mol.SetProp("smiles", Chem.MolToSmiles(mol))
                        if not mol.HasProp("inchikey"):
                            mol.SetProp("inchikey", Chem.MolToInchiKey(mol))
                        with Chem.SDWriter(f) as w:
                            w.write(mol)

    """
}

process PREPARE_LIGANDS {

    tag "${ligands}:chunk=${id}"

    publishDir(
        "${params.outputs}/ligands", 
        mode: 'copy', 
        saveAs: { v -> "${ligands.simpleName}-${id}-${v}"},
    )

    input:
    tuple val( id ), path( ligands )

    output:
    tuple val( id ), path( "mol.sdf" )

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    import sys

    from rdkit import Chem
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True  # fallback for difficult geometries

    if "${ligands}".endswith(".gz"):
        opener = partial(gzip.open, mode="rb")
    else:
        opener = partial(open, mode="rb")

    if "${ligands}".endswith((".sdf", ".sdf.gz")):
        supplier = partial(Chem.ForwardSDMolSupplier, removeHs=False)
    elif "${ligands}".endswith((".smi", ".smi.gz")):
        supplier = partial(Chem.SmilesMolSupplierFromText, titleLine=False)
    else:
        raise ValueError("Filename '${ligands}' should end with .smi or .sdf (± .gz)")

    with opener("${ligands}") as f_in:
        if "${ligands}".endswith((".sdf", ".sdf.gz")):
            supp_input = f_in
        elif "${ligands}".endswith((".smi", ".smi.gz")):
            supp_input = "\\n".join(line.decode("utf-8") for line in f_in)
        else:
            raise ValueError("Filename '${ligands}' should end with .smi or .sdf (± .gz)")
        with supplier(supp_input) as supp:
            for i, mol in enumerate(supp):
                if mol is None:
                    continue
                if mol.HasProp("zinc_id"):
                    this_name = mol.GetProp("zinc_id")
                elif mol.HasProp("name"):
                    this_name = mol.GetProp("name")
                elif mol.HasProp("_Name"):
                    this_name = mol.GetProp("_Name")
                else:
                    this_name = f"${id}-{i}"
                mol.SetProp("_Name", this_name)
                mol.SetProp("name", this_name)
                if not mol.HasProp("smiles"):
                    mol.SetProp("smiles", Chem.MolToSmiles(mol))
                if not mol.HasProp("inchikey"):
                    mol.SetProp("inchikey", Chem.MolToInchiKey(mol))
                # Standardise
                mol = Chem.AddHs(mol)
                Chem.SanitizeMol(mol)
                try:
                    result = AllChem.EmbedMolecule(mol, params)
                except Exception as e:
                    print(mol, e, file=sys.stderr)
                    result = -1
                if result == -1:
                    print(f"[WARN] mol {i}: embedding failed, skipping", file=sys.stderr)
                    continue
                AllChem.MMFFOptimizeMolecule(mol)
                with open("mol.sdf", "a") as f:
                    with Chem.SDWriter(f) as w:
                        w.write(mol)

    """
}
