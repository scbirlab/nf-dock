process DETECT_POCKET {
    /*
     * Identify the primary druggable pocket using fpocket.
     * Output: pocket center coordinates + residue list for box definition.
     *
     * If a co-crystallised ligand exists, use its centroid instead.
     */
    tag "${protein_id}"

    input:
    tuple val(protein_id), val(uniprot_id), path(structure)

    output:
    tuple val(protein_id), val(uniprot_id), path(structure), path("${protein_id}_pocket.json")

    script:
    """
    fpocket -f "${structure}"

    # Parse top-ranked pocket, extract center + residues
    python3 << 'PYEOF'
import json, glob, re
from pathlib import Path

# fpocket outputs to <name>_out/
out_dir = glob.glob("*_out")[0]
pocket_pdb = sorted(glob.glob(f"{out_dir}/pockets/pocket1_atm.pdb"))

if not pocket_pdb:
    # Fallback: whole protein centroid (blind docking)
    pocket = {"center": [0, 0, 0], "size": [30, 30, 30], "residues": []}
else:
    xs, ys, zs, residues = [], [], [], set()
    with open(pocket_pdb[0]) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                xs.append(x); ys.append(y); zs.append(z)
                resid = line[21:26].strip()
                residues.add(resid)

    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
    pad = ${params.box_padding}
    sx = max(xs) - min(xs) + 2*pad
    sy = max(ys) - min(ys) + 2*pad
    sz = max(zs) - min(zs) + 2*pad
    pocket = {
        "center": [round(cx,2), round(cy,2), round(cz,2)],
        "size":   [round(sx,2), round(sy,2), round(sz,2)],
        "residues": sorted(residues)
    }

with open("${protein_id}_pocket.json", "w") as f:
    json.dump(pocket, f)
PYEOF
    """
}