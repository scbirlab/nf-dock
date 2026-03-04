process GNINA_DOCK {
    /*
     * Core docking process. Each invocation docks one ligand chunk
     * against one protein target.
     *
     * Outputs: per-compound best Vina score, CNNscore, CNNaffinity.
     */
    tag "${protein_id}:${chunk.baseName}"
    label 'process_medium'
    cpus params.gnina_cpus
    accelerator params.gnina_gpu ? 1 : 0

    container "gnina/gnina"

    input:
    tuple val(protein_id), val(uniprot_id), path(receptor), path(pocket_json), path(chunk)

    output:
    tuple val(protein_id), val(uniprot_id), path("${protein_id}_${chunk.baseName}_scores.tsv")

    script:
    def gpu_flag = params.gnina_gpu ? "" : "--no_gpu"
    """
    # Parse pocket center and size from JSON
    CENTER=\$(python3 -c "import json; d=json.load(open('${pocket_json}')); print(f'{d[\"center\"][0]} {d[\"center\"][1]} {d[\"center\"][2]}')")
    SIZE=\$(python3 -c "import json; d=json.load(open('${pocket_json}')); print(f'{d[\"size\"][0]} {d[\"size\"][1]} {d[\"size\"][2]}')")

    CX=\$(echo \$CENTER | cut -d' ' -f1)
    CY=\$(echo \$CENTER | cut -d' ' -f2)
    CZ=\$(echo \$CENTER | cut -d' ' -f3)
    SX=\$(echo \$SIZE | cut -d' ' -f1)
    SY=\$(echo \$SIZE | cut -d' ' -f2)
    SZ=\$(echo \$SIZE | cut -d' ' -f3)

    gnina \
        -r ${receptor} \
        -l ${chunk} \
        --center_x \$CX --center_y \$CY --center_z \$CZ \
        --size_x \$SX --size_y \$SY --size_z \$SZ \
        --exhaustiveness ${params.gnina_exhaustiveness} \
        --num_modes ${params.gnina_num_modes} \
        --cnn ${params.gnina_cnn} \
        --cpu ${params.gnina_cpus} \
        ${gpu_flag} \
        -o docked.sdf.gz \
        --log scores_raw.log

    # Parse output SDF for best scores per molecule
    python3 << 'PYEOF'
import gzip
from rdkit import Chem
from collections import defaultdict

results = defaultdict(lambda: {"vina": 999, "cnnscore": 0, "cnnaffinity": 0})

suppl = Chem.ForwardSDMolSupplier(gzip.open("docked.sdf.gz"))
for mol in suppl:
    if mol is None:
        continue
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
    vina = float(mol.GetProp("minimizedAffinity")) if mol.HasProp("minimizedAffinity") else 999
    cnns = float(mol.GetProp("CNNscore")) if mol.HasProp("CNNscore") else 0
    cnna = float(mol.GetProp("CNNaffinity")) if mol.HasProp("CNNaffinity") else 0

    # Keep pose with best CNNscore (most likely to be correct)
    if cnns > results[name]["cnnscore"]:
        results[name] = {"vina": vina, "cnnscore": cnns, "cnnaffinity": cnna}

with open("${protein_id}_${chunk.baseName}_scores.tsv", "w") as f:
    f.write("compound_id\\tprotein_id\\tuniprot_id\\tvina_score\\tcnn_score\\tcnn_affinity\\n")
    for cmpd, scores in results.items():
        f.write(f"{cmpd}\\t${protein_id}\\t${uniprot_id}\\t{scores['vina']:.3f}\\t{scores['cnnscore']:.4f}\\t{scores['cnnaffinity']:.3f}\\n")
PYEOF
    """
}