"""
Part 2 -- relative frequency of each cell population per sample.

For every sample we sum the five population counts to get a total, then express
each population as a percentage of that total. The result is the long/tidy table
the assignment asks for, with exactly these columns:

    sample, total_count, population, count, percentage
"""
import pandas as pd

from analysis.db import DB_PATH, read_sql

# SQLite window function: SUM(count) OVER (PARTITION BY sample) is the per-sample
# total, computed alongside each population row in a single pass.
_FREQUENCY_SQL = """
SELECT
    sample_id                                              AS sample,
    SUM(count) OVER (PARTITION BY sample_id)               AS total_count,
    population,
    count,
    ROUND(100.0 * count / SUM(count) OVER (PARTITION BY sample_id), 4)
                                                           AS percentage
FROM cell_counts
ORDER BY sample_id, population
"""


def cell_frequencies(db_path=DB_PATH) -> pd.DataFrame:
    """Return the Part 2 relative-frequency table (one row per sample x population)."""
    return read_sql(_FREQUENCY_SQL, db_path=db_path)
