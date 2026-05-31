"""
Part 3 -- responders vs non-responders among melanoma / miraclib / PBMC samples.

We compare each population's relative frequency between responders (response
'yes') and non-responders ('no') with the Mann-Whitney U test (non-parametric:
cell fractions are bounded and not guaranteed normal). Because we run one test
per population we also report Bonferroni-adjusted p-values, and because each
subject contributes several timepoints we add a subject-level sensitivity check
that averages a subject's samples before testing.
"""
import matplotlib
matplotlib.use("Agg")  # headless: write PNGs without a display (Codespaces-safe)
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

from analysis.db import DB_PATH, POPULATIONS, read_sql

ALPHA = 0.05

# Per-sample population percentages for the cohort, with the subject metadata
# needed for grouping. Percentage is computed exactly as in Part 2.
_COHORT_SQL = """
WITH freq AS (
    SELECT
        sample_id,
        population,
        100.0 * count / SUM(count) OVER (PARTITION BY sample_id) AS percentage
    FROM cell_counts
)
SELECT
    f.sample_id        AS sample,
    f.population,
    f.percentage,
    s.subject_id,
    s.time_from_treatment_start,
    su.response
FROM freq f
JOIN samples  s  ON s.sample_id  = f.sample_id
JOIN subjects su ON su.subject_id = s.subject_id
WHERE su.condition  = 'melanoma'
  AND su.treatment  = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND su.response IN ('yes', 'no')
"""


def responder_cohort(db_path=DB_PATH) -> pd.DataFrame:
    """Long table of population percentages for the Part 3 cohort."""
    return read_sql(_COHORT_SQL, db_path=db_path)


def _mannwhitney_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Run Mann-Whitney U per population on `percentage`, grouped by response."""
    n_tests = len(POPULATIONS)
    rows = []
    for pop in POPULATIONS:
        responder = frame[(frame.population == pop) & (frame.response == "yes")].percentage
        non_resp = frame[(frame.population == pop) & (frame.response == "no")].percentage
        u_stat, p_value = stats.mannwhitneyu(responder, non_resp, alternative="two-sided")
        p_bonf = min(p_value * n_tests, 1.0)
        rows.append({
            "population": pop,
            "n_responder": len(responder),
            "n_non_responder": len(non_resp),
            "median_responder": round(responder.median(), 4),
            "median_non_responder": round(non_resp.median(), 4),
            "u_statistic": u_stat,
            "p_value": p_value,
            "p_bonferroni": p_bonf,
            "significant_p05": p_value < ALPHA,
            "significant_bonferroni": p_bonf < ALPHA,
        })
    return pd.DataFrame(rows)


def compare_responders(cohort: pd.DataFrame) -> pd.DataFrame:
    """Primary test: treat each PBMC sample as one observation."""
    return _mannwhitney_table(cohort)


def compare_responders_by_subject(cohort: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity test: average each subject's timepoints first (one row/subject)."""
    subject_means = (cohort
                     .groupby(["subject_id", "response", "population"], as_index=False)
                     .percentage.mean())
    return _mannwhitney_table(subject_means)


def significant_populations(stats_table: pd.DataFrame, bonferroni: bool = False) -> list:
    """Population names flagged significant (raw p<0.05, or Bonferroni if asked)."""
    col = "significant_bonferroni" if bonferroni else "significant_p05"
    return stats_table.loc[stats_table[col], "population"].tolist()


def boxplot(cohort: pd.DataFrame, ax=None):
    """Boxplot of relative frequency per population, responder vs non-responder."""
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(data=cohort, x="population", y="percentage", hue="response",
                order=POPULATIONS, hue_order=["yes", "no"], ax=ax)
    ax.set_title("Relative frequency: responders (yes) vs non-responders (no)\n"
                 "melanoma / miraclib / PBMC")
    ax.set_xlabel("immune cell population")
    ax.set_ylabel("relative frequency (%)")
    return ax


def save_boxplot(cohort: pd.DataFrame, path) -> None:
    """Render and write the responder/non-responder boxplot to `path`."""
    fig, ax = plt.subplots(figsize=(12, 6))
    boxplot(cohort, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
