"""Shared database access: paths, the SQLite connection, and small helpers.

Centralising this means every other module (load_data, the analyses, the
dashboard) opens the database the same way -- crucially with foreign-key
enforcement turned on, which SQLite leaves OFF by default.
"""
from pathlib import Path
import sqlite3

import pandas as pd

# Paths relative to the repo root (this file is <root>/analysis/db.py).
ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell_count.db"
OUTPUTS_DIR = ROOT / "outputs"

# The five immune-cell populations, in a fixed display order.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def connect(db_path: Path = DB_PATH) -> sqlite3.Connection:
    """Open the SQLite database with foreign-key constraints enforced."""
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def read_sql(query: str, params=None, db_path: Path = DB_PATH) -> pd.DataFrame:
    """Run a SELECT and return a DataFrame (connection opened/closed for you)."""
    with connect(db_path) as conn:
        return pd.read_sql_query(query, conn, params=params)
