process PocketDetect_P2Rank {

    tag "${id}-${uniprot_id}"

    publishDir(
        "${params.outputs}/pockets", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${v}"},
    )

    input:
    tuple val( id ), val( uniprot_id ), path( structure )

    output:
    tuple val( id ), val( uniprot_id ), path( structure ), path( "pocket.json" ), emit: main
    tuple val( id ), val( uniprot_id ), path( "p2rank_out/*_predictions.csv" ), emit: info

    script:
    """
    prank predict -f "${structure}" -o p2rank_out -c alphafold
    mv p2rank_out/*.csv info.csv

    # Parse top-ranked pocket, extract center + residues
    python -c '
    import csv
    from collections import defaultdict
    from glob import glob
    import gzip
    import json
    import os
    import numpy as np

    base = os.path.basename("${structure}")
    # P2Rank appends _predictions.csv to the full filename
    pred_csv  = f"p2rank_out/{base}_predictions.csv"
    points_gz = f"p2rank_out/visualizations/data/{base}_points.pdb.gz"

    # --- Parse predictions CSV ---
    pockets_meta = []
    if os.path.exists(pred_csv):
        with open(pred_csv) as f:
            reader = csv.DictReader(f, skipinitialspace=True)
            for row in reader:
                pockets_meta.append(row)
    pockets_meta = sorted(pockets_meta, key=lambda r: int(r["rank"]))[:2]

    if not pockets_meta:
        # Fallback: whole protein centroid (blind docking)
        pocket_list = [
            {"center": [0, 0, 0], "size": [30, 30, 30], "residues": [], "blind": True, "id": 0}
        ]
    else:
        rank_coords = defaultdict(list)
        if os.path.exists(points_gz):
            with gzip.open(points_gz, "rt") as f:
                for line in f:
                    rank = int(line[22:26].strip())
                    if line.startswith(("ATOM","HETATM")) and rank > 0:
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        rank_coords[rank].append([x, y, z])
        pocket_list = []
        for meta in pockets_meta:
            rank = int(meta["rank"])
            cx = float(meta["center_x"])
            cy = float(meta["center_y"])
            cz = float(meta["center_z"])
            if rank in rank_coords:
                pad = ${params.box_padding}
                pts = np.asarray(rank_coords[rank])
                sx, sy, sz = (
                    max(20., float(pts[:,i].max() - pts[:,i].min()) + 2. * pad)
                    for i in range(3)
                )
            else:
                sx = sy = sz = 20.
            residues = [r.strip() for r in meta.get("residue_ids","").split() if r.strip()]
            pocket_list.append({
                "center": [cx, cy, cz],
                "size": [sx, sy, sz],
                "residues": sorted(residues),
                "blind": False,
                "id": rank - 1,
            })

    with open("pocket.json", "w") as f:
        json.dump(pocket_list, f, indent=4)
    '

    """
}
