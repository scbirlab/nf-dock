process PREPARE_LIGANDS {
    /*
     * Standardise, generate 3D conformers, protonate at pH 7.4.
     * Split into individual mol2/SDF chunks for parallel docking.
     */

    input:
    path(ligands_sdf)

    output:
    path("chunks/*.sdf")

    script:
    """
    mkdir -p chunks

    python3 << 'PYEOF'
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import os, math

suppl = Chem.SDMolSupplier("${ligands_sdf}", removeHs=False)
mols = []
for mol in suppl:
    if mol is None:
        continue
    # Standardise
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    mols.append(mol)

# Chunk into batches of 100 for parallel docking
chunk_size = 100
n_chunks = math.ceil(len(mols) / chunk_size)
for i in range(n_chunks):
    chunk = mols[i*chunk_size : (i+1)*chunk_size]
    writer = Chem.SDWriter(f"chunks/chunk_{i:05d}.sdf")
    for m in chunk:
        writer.write(m)
    writer.close()
PYEOF
    """
}
