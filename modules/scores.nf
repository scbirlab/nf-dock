process AGGREGATE_SCORES {
    /*
     * Concatenate all per-chunk score files into one matrix.
     * Output: compounds × proteins score matrix (long format + wide pivot).
     */
    publishDir "${params.outdir}/docking", mode: 'copy'

    input:
    path(score_files)    // collected list of all TSV files

    output:
    path("score_matrix_long.parquet"), emit: scores_long
    path("score_matrix_wide.parquet"), emit: scores_wide

    script:
    """
    python3 << 'PYEOF'
import pandas as pd
import glob

files = glob.glob("*_scores.tsv")
dfs = [pd.read_csv(f, sep="\\t") for f in files]
long = pd.concat(dfs, ignore_index=True)
long.to_parquet("score_matrix_long.parquet")

# Pivot: compounds as rows, proteins as columns, CNNaffinity as values
wide = long.pivot_table(
    index="compound_id",
    columns="protein_id",
    values="cnn_affinity",
    aggfunc="max"
)
wide.to_parquet("score_matrix_wide.parquet")
PYEOF
    """
}