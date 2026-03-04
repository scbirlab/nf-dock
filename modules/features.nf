process COMPUTE_FEATURES {
    /*
     * Pre-compute compound fingerprints and protein embeddings
     * for surrogate model training.
     *
     * Compound: Morgan fingerprints (radius 2, 2048 bits)
     * Protein:  ESM-2 mean-pooled embeddings (cached)
     */
    conda "conda-forge::rdkit pip::fair-esm"
    label 'process_high'
    accelerator 1   // GPU for ESM-2

    input:
    path(ligands_sdf)
    path(protein_fastas)   // FASTA of all protein sequences

    output:
    path("compound_fps.npz"),   emit: compound_fps
    path("protein_embeds.npz"), emit: protein_embeds

    script:
    """
    python3 << 'PYEOF'
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# ── Compound fingerprints ──
suppl = Chem.SDMolSupplier("${ligands_sdf}", removeHs=True)
names, fps = [], []
for mol in suppl:
    if mol is None: continue
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{len(names)}"
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=${params.fp_radius}, nBits=${params.fp_nbits})
    names.append(name)
    fps.append(np.array(fp))

np.savez("compound_fps.npz", names=names, fps=np.array(fps))

# ── Protein embeddings via ESM-2 ──
import torch, esm
from Bio import SeqIO

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval().cuda()
batch_converter = alphabet.get_batch_converter()

protein_names, protein_embeds = [], []
for record in SeqIO.parse("${protein_fastas}", "fasta"):
    seq = str(record.seq)[:1022]  # ESM-2 max length
    data = [(record.id, seq)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.cuda()

    with torch.no_grad():
        result = model(tokens, repr_layers=[33])
    # Mean pool over sequence (exclude BOS/EOS)
    embed = result["representations"][33][0, 1:len(seq)+1].mean(dim=0).cpu().numpy()
    protein_names.append(record.id)
    protein_embeds.append(embed)

np.savez("protein_embeds.npz", names=protein_names, embeds=np.array(protein_embeds))
PYEOF
    """
}