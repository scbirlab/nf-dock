

process TRAIN_SURROGATE {
    /*
     * Train surrogate model: f(compound_fp, protein_embed) → CNNaffinity
     *
     * Uses the docked subset (Phase 2 output) as training data.
     * Evaluates on held-out compound-protein pairs.
     */
    publishDir "${params.outdir}/model", mode: 'copy'
    conda "conda-forge::scikit-learn conda-forge::xgboost conda-forge::pandas"

    input:
    path(scores_long)
    path(compound_fps)
    path(protein_embeds)

    output:
    path("surrogate_model.pkl"),    emit: model
    path("training_metrics.json"),  emit: metrics

    script:
    """
    python3 << 'PYEOF'
import numpy as np, pandas as pd, json, pickle
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_squared_error, r2_score

# Load
scores = pd.read_parquet("${scores_long}")
cfp = np.load("${compound_fps}", allow_pickle=True)
pem = np.load("${protein_embeds}", allow_pickle=True)

cpd_idx = {n: i for i, n in enumerate(cfp["names"])}
prt_idx = {n: i for i, n in enumerate(pem["names"])}

# Build feature matrix: [compound_fp | protein_embed]
X, y, groups = [], [], []
for _, row in scores.iterrows():
    ci = cpd_idx.get(row["compound_id"])
    pi = prt_idx.get(row["protein_id"])
    if ci is None or pi is None:
        continue
    feat = np.concatenate([cfp["fps"][ci], pem["embeds"][pi]])
    X.append(feat)
    y.append(row["cnn_affinity"])
    groups.append(row["compound_id"])  # group-split by compound

X, y = np.array(X), np.array(y)

# Train with grouped CV (no compound leakage)
if "${params.surrogate_type}" == "xgboost":
    from xgboost import XGBRegressor
    model = XGBRegressor(
        n_estimators=500, max_depth=8, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        tree_method="hist", n_jobs=-1
    )
else:
    from sklearn.neural_network import MLPRegressor
    model = MLPRegressor(
        hidden_layer_sizes=(512, 256, 128), max_iter=200,
        early_stopping=True, validation_fraction=0.1
    )

# Simple train/test split by compound groups
gkf = GroupKFold(n_splits=5)
r2s, rmses = [], []
for train_idx, test_idx in gkf.split(X, y, groups):
    model.fit(X[train_idx], y[train_idx])
    pred = model.predict(X[test_idx])
    r2s.append(r2_score(y[test_idx], pred))
    rmses.append(mean_squared_error(y[test_idx], pred, squared=False))

# Final model on all data
model.fit(X, y)
with open("surrogate_model.pkl", "wb") as f:
    pickle.dump(model, f)

metrics = {
    "cv_r2_mean": float(np.mean(r2s)),
    "cv_r2_std":  float(np.std(r2s)),
    "cv_rmse_mean": float(np.mean(rmses)),
    "n_train": len(y)
}
with open("training_metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)

print(f"CV R²: {np.mean(r2s):.3f} ± {np.std(r2s):.3f}")
print(f"CV RMSE: {np.mean(rmses):.3f}")
PYEOF
    """
}


process PREDICT_FULL_MATRIX {
    /*
     * Use surrogate to predict docking scores for all
     * compound-protein pairs (including those not explicitly docked).
     */
    publishDir "${params.outdir}/predictions", mode: 'copy'
    conda "conda-forge::scikit-learn conda-forge::xgboost conda-forge::pandas"

    input:
    path(model)
    path(compound_fps)
    path(protein_embeds)

    output:
    path("full_matrix.parquet"), emit: full_matrix

    script:
    """
    python3 << 'PYEOF'
import numpy as np, pandas as pd, pickle

with open("${model}", "rb") as f:
    model = pickle.load(f)

cfp = np.load("${compound_fps}", allow_pickle=True)
pem = np.load("${protein_embeds}", allow_pickle=True)

# Predict all pairs
rows = []
for pi, pname in enumerate(pem["names"]):
    # Batch: tile protein embed across all compounds
    P = np.tile(pem["embeds"][pi], (len(cfp["fps"]), 1))
    X = np.concatenate([cfp["fps"], P], axis=1)
    preds = model.predict(X)

    for ci, score in enumerate(preds):
        rows.append({
            "compound_id": cfp["names"][ci],
            "protein_id": pname,
            "predicted_affinity": float(score)
        })

df = pd.DataFrame(rows)
df.to_parquet("full_matrix.parquet")
PYEOF
    """
}
