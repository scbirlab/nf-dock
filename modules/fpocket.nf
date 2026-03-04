process DETECT_POCKET {
    /*
     * Identify the primary druggable pocket using fpocket.
     * Output: pocket center coordinates + residue list for box definition.
     *
     * If a co-crystallised ligand exists, use its centroid instead.
     */
    tag "${id}"

    publishDir(
        "${params.outputs}/pockets", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${v}"},
    )

    input:
    tuple val( id ), val( uniprot_id ), path( structure )

    output:
    tuple val( id ), val( uniprot_id ), path( structure ), path( "pocket.json" ), emit: main
    tuple val( id ), val( uniprot_id ), path( "info.txt" ), emit: info

    script:
    """
    fpocket -f "${structure}"
    mv *_out/*_info.txt info.txt

    # Parse top-ranked pocket, extract center + residues
    python -c '
    from glob import glob
    import json
    import os

    # fpocket outputs to <name>_out/
    out_dir = "*_out"
    pocket_pdb = sorted(glob(os.path.join(out_dir, "pockets", "pocket*_atm.pdb")))[:3]
    #pocket_pdb = sorted(glob(os.path.join(out_dir, "*_out.pdb")))

    pocket_list = []
    if not pocket_pdb:
        # Fallback: whole protein centroid (blind docking)
        pocket_list = [
            {"center": [0, 0, 0], "size": [30, 30, 30], "residues": [], "blind": True, "id": 0}
        ]
    else:
        xs, ys, zs, residues = [], [], [], set()
        for i, pocket in enumerate(pocket_pdb):
            with open(pocket) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        xs.append(x); ys.append(y); zs.append(z)
                        resid = line[21:26].strip()
                        residues.add(resid)

            cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
            pad = ${params.box_padding}
            sx, sy, sz = (
                max(dim) - min(dim) + 2. * pad
                for dim in (xs, ys, zs)
            )
            pocket_list.append({
                "center": [cx, cy, cz],
                "size":   [sx, sy, sz],
                "residues": sorted(residues),
                "blind": False,
                "id": i,
            })

    with open("pocket.json", "w") as f:
        json.dump(pocket_list, f, indent=4)
    '

    """
}

process SplitPockets {

    tag "${id}"

    input:
    tuple val( id ), val( uniprot_id ), path( structure ), path( pocket )

    output:
    tuple val( id ), val( uniprot_id ), path( structure ), path( "pocket_*.json" )

    script:
    """
    #!/usr/bin/env python

    import json

    with open("${pocket}", "r") as f:
        d = json.load(f)
        for i, pocket in enumerate(d):
            with open(f"pocket_{i:02d}.json", "w") as o:
                json.dump(pocket, o, indent=4)
    
    """

}