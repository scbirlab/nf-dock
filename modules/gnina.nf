process GNINA_DOCK {
    
    tag "${id}:${pocket_json.baseName}"
    label 'gpu_single'

    publishDir(
        "${params.outputs}/docking", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${pocket_json.baseName}-${v}"},
    )

    container "gnina/gnina"

    input:
    tuple val( id ), val( uniprot_id ), path(receptor), path( pocket_json ), path( chunk )

    output:
    tuple val( id ), val( uniprot_id ), path( pocket_json ), path( "docked.sdf*" ), emit: main
    tuple val( id ), path( "gnina.log" ), emit:logs

    script:
    """
    set -euox pipefail
    # Parse pocket center and size from JSON
    CENTER=\$(python3 -c 'import json; d = json.load(open("${pocket_json}")); print(*d["center"])')
    SIZE=\$(python3 -c 'import json; d = json.load(open("${pocket_json}")); print(*d["size"])')

    CX=\$(echo \$CENTER | cut -d' ' -f1)
    CY=\$(echo \$CENTER | cut -d' ' -f2)
    CZ=\$(echo \$CENTER | cut -d' ' -f3)
    SX=\$(echo \$SIZE | cut -d' ' -f1)
    SY=\$(echo \$SIZE | cut -d' ' -f2)
    SZ=\$(echo \$SIZE | cut -d' ' -f3)

    cat ${chunk} > concat_chunks.sdf
    #obabel concat_chunks.sdf -O concat_cleaned.sdf

    gnina \
        -r "${receptor}" \
        -l concat_chunks.sdf \
        --center_x \$CX --center_y \$CY --center_z \$CZ \
        --size_x \$SX --size_y \$SY --size_z \$SZ \
        --exhaustiveness "${params.gnina_exhaustiveness}" \
        --num_modes "${params.gnina_num_modes}" \
        --scoring vinardo \
        --cnn_scoring rescore \
        --cnn "${params.gnina_cnn}" \
        --pose_sort_order CNNaffinity \
        --atom_term_data \
        --cpu "${task.cpus}" \
        -o docked.sdf \
        --seed 42 \
        --log gnina.log
    
    """
}

process Extract_Gnina_scores {

    tag "${id}:${chunk_id}:${pocket_json.baseName}"

    input:
    tuple val( id ), val( uniprot_id ), path(pocket_json), path( docking )

    output:
    tuple val( id ), path( "scores.tsv" )

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    import gzip
    import json

    import pandas as pd
    from rdkit import Chem

    
    with open("${pocket_json}") as f:
        d = json.load(f)

    opener = partial(gzip.open if "${docking}".endswith(".gz") else open, mode="rb")

    results = []
    with opener("${docking}") as f:
        for i, mol in enumerate(Chem.ForwardSDMolSupplier(f)):
            if mol is None:
                continue
            prop_dict = mol.GetPropsAsDict()
            results.append({
                "ligand_name": prop_dict.get("_Name"),
                "ligand_smiles": prop_dict.get("smiles"),
                "ligand_zinc_id": prop_dict.get("zinc_id"),
                "receptor_uniprot_id": "${uniprot_id}",
                "receptor_pocket_id": d["id"],
                "receptor_pocket_center": d["center"],
                "receptor_pocket_size": d["size"],
                "blind_dock": d["blind"],
                "pose_id": f"{i:06d}", 
                "affinity": float(prop_dict.get("minimizedAffinity")), 
                "cnn_score": float(prop_dict.get("CNNscore")), 
                "cnn_affinity": float(prop_dict.get("CNNaffinity")),
                "cnn_affinity_var": float(prop_dict.get("CNNaffinity_variance", 0.)),
                "cnn_vs": float(prop_dict.get("CNN_VS")),
            } | prop_dict)

    pd.DataFrame(
        results,
    ).to_csv("scores.tsv", index=False, sep="\\t")

    """
}

process AGGREGATE_SCORES {
    /*
     * Concatenate all per-chunk score files into one matrix.
     * Output: compounds × proteins score matrix (long format + wide pivot).
     */
    publishDir "${params.outputs}/docking", mode: 'copy'

    input:
    path ( score_files, stageAs: 'score_????????/scores.tsv' )

    output:
    path "score_matrix_long.tsv", emit: scores_long
    path "score_matrix_wide.tsv", emit: scores_wide

    script:
    """
    #!/usr/bin/env python

    from glob import glob
    import os

    import pandas as pd

    files = glob(os.path.join("score_*", "scores.tsv"))
    dfs = [pd.read_csv(f, sep="\\t") for f in files]
    long = pd.concat(dfs, axis=0)
    long.to_csv("score_matrix_long.tsv", sep="\\t", index=False)

    # Pivot: compounds as rows, proteins as columns, CNNaffinity as values
    ligand_cols = ["ligand_name", "ligand_smiles", "ligand_zinc_id"]
    receptor_cols = ["receptor_uniprot_id"]
    value = "cnn_affinity"
    wide = (
        long
        .sort_values(value)
        .groupby(ligand_cols + receptor_cols, dropna=False)
        .tail(1)
        .pivot(
            index=ligand_cols,
            columns=receptor_cols[0],
            values=value,
        )
    )
    wide.to_csv("score_matrix_wide.tsv", sep="\\t")

    """
}