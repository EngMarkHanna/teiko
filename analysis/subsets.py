"""
Part 4 -- baseline subset analysis.

The subset of interest: melanoma PBMC samples at baseline
(time_from_treatment_start = 0) from subjects treated with miraclib. From that
subset we break down samples by project, subjects by response and by sex, and
compute the average B-cell count for male responders.
"""
import pandas as pd

from analysis.db import DB_PATH, read_sql

# One row per baseline sample, carrying the subject metadata and the b_cell
# count we need (pulled from the long cell_counts table via the population key).
_BASELINE_SQL = """
SELECT
    s.sample_id,
    su.subject_id,
    su.project_id,
    su.response,
    su.sex,
    bc.count AS b_cell
FROM samples s
JOIN subjects su   ON su.subject_id = s.subject_id
JOIN cell_counts bc ON bc.sample_id = s.sample_id AND bc.population = 'b_cell'
WHERE su.condition  = 'melanoma'
  AND su.treatment  = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND s.time_from_treatment_start = 0
"""


def baseline_subset(db_path=DB_PATH) -> pd.DataFrame:
    """Melanoma / miraclib / PBMC / baseline samples with b_cell counts."""
    return read_sql(_BASELINE_SQL, db_path=db_path)


def samples_per_project(subset: pd.DataFrame) -> pd.DataFrame:
    """How many samples come from each project."""
    return (subset.groupby("project_id").size()
            .reset_index(name="n_samples").sort_values("project_id"))


def subjects_by_response(subset: pd.DataFrame) -> pd.DataFrame:
    """How many subjects were responders vs non-responders."""
    subjects = subset.drop_duplicates("subject_id")
    return (subjects.groupby("response").size()
            .reset_index(name="n_subjects").sort_values("response"))


def subjects_by_sex(subset: pd.DataFrame) -> pd.DataFrame:
    """How many subjects were male vs female."""
    subjects = subset.drop_duplicates("subject_id")
    return (subjects.groupby("sex").size()
            .reset_index(name="n_subjects").sort_values("sex"))


def male_responder_mean_bcells(subset: pd.DataFrame) -> float:
    """Average B-cell count for male responders at baseline (two decimals)."""
    males_resp = subset[(subset.sex == "M") & (subset.response == "yes")]
    return round(males_resp.b_cell.mean(), 2)
