process GNINA_DOCK {
    
    tag "${id}:${pocket_json.baseName}:chunk=${chunk_id}"
    label 'big_cpu'
    time '7d'

    errorStrategy "ignore" //{ ( task.attempt > 2 ) ? 'ignore' : 'retry' }

    publishDir(
        "${params.outputs}/docking/logs", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${pocket_json.baseName}-${chunk_id}-${v}" },
        pattern: "*.log"
    )

    container "gnina/gnina"

    input:
    tuple val( id ), val( uniprot_id ), path( receptor ), path( pocket_json ), val( chunk_id ), path( chunk )

    output:
    tuple val( id ), val( uniprot_id ), val( chunk_id ), path( pocket_json ), path( "docked.sdf*" ), emit: main
    tuple val( id ), val( chunk_id ), path( "gnina.log" ), emit:logs

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

    gnina \
        -r "${receptor}" \
        -l "${chunk}" \
        --center_x \$CX --center_y \$CY --center_z \$CZ \
        --size_x \$SX --size_y \$SY --size_z \$SZ \
        --exhaustiveness "${params.gnina_exhaustiveness}" \
        --num_modes "${params.gnina_num_modes}" ${params.gnina_minimize ? "--minimize" : ""} \
        --scoring vinardo \
        --cnn_scoring rescore \
        --cnn "${params.gnina_cnn}" \
        --pose_sort_order CNNaffinity \
        --cpu "${task.cpus}" \
        -o docked.sdf \
        --seed 42 \
        --log gnina.log

    gzip --best docked.sdf
    
    """
}


process Concatenate_ligands {

    tag "${id}:${pocket_json.baseName}"

    container "gnina/gnina"

    input:
    tuple val( id ), val( uniprot_id ), path( pocket_json ), path( sdfs )

    output:
    tuple val( id ), val( uniprot_id ), path( pocket_json ), path( "all.sdf.gz" )

    script:
    """
    set -euox pipefail

    zcat ${sdfs} | gzip --best > all.sdf.gz

    """
}

process Extract_ligand_receptor_poses {

    tag "${id}:${pocket_json.baseName}"

    publishDir(
        "${params.outputs}/docking/ligands", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}/${pocket_json.baseName}-${v}" },
        pattern: "*.{sdf,log}"
    )

    input:
    tuple val( id ), val( uniprot_id ), path( pocket_json ), path( ligands )

    output:
    tuple val( id ), val( uniprot_id ), path( pocket_json ), path( "*-*-*.sdf*" )

    script:
    """
    #!/usr/bin/env python

    from functools import partial
    import gzip
    from itertools import batched

    from rdkit import Chem

    supplier = partial(Chem.ForwardSDMolSupplier, removeHs=False)
    if "${ligands}".endswith(".gz"):
        opener = partial(gzip.open, mode="rb")
    else:
        opener = partial(open, mode="rb")

    with supplier(opener("${ligands}")) as supp:
        for i, mol in enumerate(supp):
            name = mol.GetProp("name")
            inchikey = mol.GetProp("inchikey")
            with open(f"{name}-{inchikey}.sdf", "a") as f:
                with Chem.SDWriter(f) as w:
                    w.write(mol)
    
    """


}

process Extract_Gnina_scores {

    tag "${id}:${pocket_json.baseName}:chunk=${chunk_id}"

    publishDir(
        "${params.outputs}/docking/scores", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${pocket_json.baseName}-${chunk_id}-${v}" },
        pattern: "*.{sdf,log}"
    )

    input:
    tuple val( id ), val( uniprot_id ), val( chunk_id ), path( pocket_json ), path( docking )

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
                "ligand_name": prop_dict.get("name"),
                "ligand_smiles": prop_dict.get("smiles"),
                "ligand_inchikey": prop_dict.get("inchikey"),
                "ligand_zinc_id": prop_dict.get("zinc_id"),
                "receptor_id": "${id}",
                "receptor_uniprot_id": "${uniprot_id}",
                "receptor_pocket_id": d["id"],
                "receptor_pocket_center": d["center"],
                "receptor_pocket_size": d["size"],
                "blind_dock": d["blind"],
                "pose_id": f"{i:06d}",
                "library_chunk_id": int(${chunk_id}),
                "affinity": float(prop_dict.get("minimizedAffinity")), 
                "cnn_score": float(prop_dict.get("CNNscore")), 
                "cnn_affinity": float(prop_dict.get("CNNaffinity")),
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
    path "score_matrix_wide_z*.tsv", emit: scores_wide_z

    script:
    """
    #!/usr/bin/env python

    from glob import glob
    import os

    import numpy as np
    import pandas as pd
    from scipy.stats import median_abs_deviation

    def zscore(m: pd.DataFrame, axis=0):
        ma = np.asarray(m)
        m_mean = np.median(ma, axis=axis, keepdims=True)
        m_std = median_abs_deviation(ma, axis=axis, keepdims=True)
        return (m - m_mean) / m_std

    files = glob(os.path.join("score_*", "scores.tsv"))
    dfs = [pd.read_csv(f, sep="\\t") for f in files]
    long = pd.concat(dfs, axis=0)
    long.to_csv("score_matrix_long.tsv", sep="\\t", index=False)

    # Pivot: compounds as rows, proteins as columns, CNNaffinity as values
    ligand_cols = ["ligand_name", "ligand_smiles", "ligand_inchikey", "ligand_zinc_id"]
    receptor_cols = ["receptor_id", "receptor_uniprot_id"]
    value = "cnn_affinity"
    wide = (
        long
        .sort_values(value)
        .groupby(ligand_cols + receptor_cols, dropna=False)
        .tail(1)
        .assign(
            receptor_id=lambda x: x["receptor_id"].str.cat(x["receptor_uniprot_id"], sep=":"),
        )
        .pivot(
            index=ligand_cols,
            columns="receptor_id",
            values=value,
        )
    )
    wide.to_csv("score_matrix_wide.tsv", sep="\\t")

    wide_z0 = zscore(wide, axis=1)
    wide_z0.to_csv("score_matrix_wide_z-col.tsv", sep="\\t")
    wide_z = zscore(wide_z0, axis=0)
    wide_z.to_csv("score_matrix_wide_z.tsv", sep="\\t")

    """
}
