process ASSIGN_TARGETS {
    /*
     * For each compound, identify the best-scoring bacterial protein
     * ("target assignment"). Then use ortholog table to predict
     * which species would be inhibited.
     */
    publishDir "${params.outdir}/atlas", mode: 'copy'

    input:
    path(full_matrix)
    path(ortholog_table)

    output:
    path("target_assignments.tsv"),  emit: assignments
    path("selectivity_matrix.tsv"),  emit: selectivity

    script:
    """
    python3 << 'PYEOF'
import pandas as pd
import numpy as np

# ── Target assignment ──
matrix = pd.read_parquet("${full_matrix}")

# Best target = protein with highest predicted affinity for each compound
best = matrix.loc[matrix.groupby("compound_id")["predicted_affinity"].idxmax()]
best = best.rename(columns={"protein_id": "assigned_target", "predicted_affinity": "best_score"})

# Confidence: gap between best and second-best score
def score_gap(grp):
    top2 = grp.nlargest(2, "predicted_affinity")["predicted_affinity"].values
    return top2[0] - top2[1] if len(top2) > 1 else np.nan

gaps = matrix.groupby("compound_id").apply(score_gap).rename("score_gap")
best = best.merge(gaps, on="compound_id")
best.to_csv("target_assignments.tsv", sep="\\t", index=False)

# ── Selectivity prediction via orthology ──
ortho = pd.read_csv("${ortholog_table}", sep="\\t")
# Columns: ref_protein_id, species, ortholog_protein_id, aa_identity

# For each compound's assigned target, check which species have orthologs
merged = best.merge(ortho, left_on="assigned_target", right_on="ref_protein_id", how="left")

# Predict: species with high-identity ortholog → likely susceptible
# Threshold could be tuned; using 40% as in RepOrt
merged["predicted_susceptible"] = merged["aa_identity"] >= 40.0

selectivity = merged.pivot_table(
    index="compound_id",
    columns="species",
    values="predicted_susceptible",
    aggfunc="max",
    fill_value=False
)
selectivity.to_csv("selectivity_matrix.tsv", sep="\\t")
PYEOF
    """
}

process VALIDATE {
    /*
     * Compare predictions against known antibacterial activities
     * (e.g., Maier 2018 dataset).
     */
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path(selectivity)
    path(known_activities)   // TSV: compound_id, species, active (bool)

    output:
    path("validation_report.json"), emit: report

    script:
    """
    python3 << 'PYEOF'
import pandas as pd, json
from sklearn.metrics import roc_auc_score, precision_score, recall_score

pred = pd.read_csv("${selectivity}", sep="\\t", index_col=0)
known = pd.read_csv("${known_activities}", sep="\\t")

# Align
results = []
for _, row in known.iterrows():
    cpd, sp, active = row["compound_id"], row["species"], row["active"]
    if cpd in pred.index and sp in pred.columns:
        results.append({
            "compound_id": cpd, "species": sp,
            "actual": active, "predicted": pred.loc[cpd, sp]
        })

rdf = pd.DataFrame(results)
if len(rdf) > 0 and rdf["actual"].nunique() > 1:
    metrics = {
        "n_pairs": len(rdf),
        "auroc": float(roc_auc_score(rdf["actual"], rdf["predicted"])),
        "precision": float(precision_score(rdf["actual"], rdf["predicted"])),
        "recall": float(recall_score(rdf["actual"], rdf["predicted"]))
    }
else:
    metrics = {"n_pairs": len(rdf), "note": "insufficient data for metrics"}

with open("validation_report.json", "w") as f:
    json.dump(metrics, f, indent=2)

print(json.dumps(metrics, indent=2))
PYEOF
    """
}
