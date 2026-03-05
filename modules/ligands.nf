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
    if ( ligands.extension == "smi" ) {
        """
        split -l${batch_size} -d "${ligands}" lib_ --additional-suffix=".smi"

        """
    }
    else {
        """
        #!/usr/bin/env python

        from itertools import batched

        from rdkit import Chem
        from rdkit.Chem import AllChem

        with Chem.SDMolSupplier("${ligands}") as supp:
            for i, mols in enumerate(batched(supp, ${batch_size})):
                with open(f"lib_{i:06d}.sdf", "a") as f:
                    for j, mol in enumerate(mols):
                        if mol is None:
                            continue
                        if not mol.HasProp("_Name"):
                            mol.SetProp("_Name", f"{i}-{j}")
                        mol.SetProp("smiles", Chem.MolToSmiles(mol))
                        with Chem.SDWriter(f) as w:
                            w.write(mol)

        """
    }
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

    if "${ligands}".endswith((".sdf", ".sdf.gz")):
        opener = partial(Chem.SDMolSupplier, removeHs=False)
    elif "${ligands}".endswith((".smi", ".smi.gz")):
        opener = partial(Chem.SmilesMolSupplier, titleLine=False)

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True  # fallback for difficult geometries

    with opener("${ligands}") as supp:
        for i, mol in enumerate(supp):
            if mol is None:
                continue
            if not mol.HasProp("_Name"):
                mol.SetProp("_Name", f"${id}-{i}")
            if not mol.HasProp("smiles"):
                mol.SetProp("smiles", Chem.MolToSmiles(mol))
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
