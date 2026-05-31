"""
load_data.py  --  Part 1: initialise the SQLite database and load the CSV.

Run from the repository root, no arguments:

    python load_data.py

It creates cell_count.db, applies schema.sql, and loads every row of
cell-count.csv into the four normalised tables (projects, subjects, samples,
cell_counts). Re-running is safe: schema.sql drops and recreates the tables.
"""
import sys

import pandas as pd

from analysis.db import (CSV_PATH, DB_PATH, POPULATIONS, SCHEMA_PATH, connect)


def build_frames(df: pd.DataFrame):
    """Split one flat CSV into the four table-shaped DataFrames."""
    projects = (df[["project"]]
                .drop_duplicates()
                .rename(columns={"project": "project_id"}))

    subjects = (df[["subject", "project", "condition", "age",
                    "sex", "treatment", "response"]]
                .drop_duplicates(subset="subject")
                .rename(columns={"subject": "subject_id",
                                 "project": "project_id"}))

    samples = (df[["sample", "subject", "sample_type",
                   "time_from_treatment_start"]]
               .rename(columns={"sample": "sample_id",
                                "subject": "subject_id"}))

    # Wide -> long: 5 population columns become (sample_id, population, count).
    cell_counts = (df.melt(id_vars="sample", value_vars=POPULATIONS,
                           var_name="population", value_name="count")
                   .rename(columns={"sample": "sample_id"}))

    return projects, subjects, samples, cell_counts


def main() -> int:
    print(f"Reading {CSV_PATH.name} ...")
    df = pd.read_csv(CSV_PATH)
    print(f"  {len(df):,} rows")

    projects, subjects, samples, cell_counts = build_frames(df)

    print(f"Creating database {DB_PATH.name} from {SCHEMA_PATH.name} ...")
    conn = connect()
    try:
        conn.executescript(SCHEMA_PATH.read_text())

        # Insert parents before children so foreign keys always resolve.
        for name, frame in [("projects", projects), ("subjects", subjects),
                            ("samples", samples), ("cell_counts", cell_counts)]:
            frame.to_sql(name, conn, if_exists="append", index=False)
            print(f"  loaded {name:<12} {len(frame):>7,} rows")
        conn.commit()

        # Confirm what actually landed in the database.
        print("Verifying row counts in database ...")
        ok = True
        expected = {"projects": len(projects), "subjects": len(subjects),
                    "samples": len(samples), "cell_counts": len(cell_counts)}
        for name, exp in expected.items():
            got = conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
            status = "ok" if got == exp else "MISMATCH"
            ok = ok and got == exp
            print(f"  {name:<12} {got:>7,}  ({status})")
        if not ok:
            print("ERROR: row-count mismatch after load", file=sys.stderr)
            return 1
    finally:
        conn.close()

    print(f"\nDone. Database ready at {DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
