# Loblaw Bio — Immune Cell Population Analysis

A small, reproducible pipeline that loads the clinical-trial cell-count data into
a normalised SQLite database, answers the four analysis questions, and serves the
results in an interactive dashboard.

- **Part 1** — relational schema + loader (`schema.sql`, `load_data.py`)
- **Part 2** — relative frequency of each population per sample
- **Part 3** — responders vs non-responders (melanoma / miraclib / PBMC)
- **Part 4** — baseline subset breakdowns
- **Dashboard** — Streamlit app mirroring Parts 2–4

---

## Quick start (GitHub Codespaces or local)

```bash
make setup       # create .venv and install dependencies
make pipeline    # build the database and generate all tables + plots
make dashboard   # launch the interactive dashboard
```

`make pipeline` runs the whole thing with no manual steps: it builds
`cell_count.db` from `cell-count.csv`, then computes Parts 2–4, writes everything
to `outputs/`, and runs regression checks that fail loudly if any verified number
drifts.

> Requires Python 3.11+ and `make`. SQLite needs **no installation** — it is part
> of the Python standard library (`sqlite3`); the database is simply the file
> `cell_count.db` created on first connect.

**Dashboard link:** [local/Codespaces dashboard](http://localhost:8501) after
running `make dashboard`. In Codespaces, open the forwarded port 8501.

For the final public submission link, deploy this repo on
[Streamlit Community Cloud](https://streamlit.io/cloud) with `dashboard.py` as the
entrypoint, then replace the local link above with the deployed app URL.

---

## Repository layout

```
.
├── cell-count.csv              input data (given)
├── schema.sql                  DDL: the 4 tables + indexes (single source of truth)
├── load_data.py                Part 1: create cell_count.db and load every row
├── run_pipeline.py             runs Parts 2–4, writes outputs/, runs regression checks
├── analysis/                   reusable analysis logic (imported by pipeline & dashboard)
│   ├── db.py                   paths + SQLite connection (foreign keys ON)
│   ├── frequencies.py          Part 2: relative-frequency table
│   ├── stats_tests.py          Part 3: cohort, Mann-Whitney, boxplot
│   └── subsets.py              Part 4: baseline subset queries
├── dashboard.py                Streamlit dashboard (Parts 2–4)
├── outputs/                    generated tables (.csv) and plots (.png)
├── datanalysis/                exploratory analysis used to design the schema
├── requirements.txt            pinned dependencies
├── Makefile                    setup / pipeline / dashboard
└── README.md
```

### Why this structure

- **`load_data.py` is in the root and takes no arguments** — run it directly with
  `python load_data.py`, exactly as required.
- **One module per analysis Part**, collected in a small `analysis/` package. Each
  function returns a tidy DataFrame and is imported by **both** `run_pipeline.py`
  (batch outputs) and `dashboard.py` (interactive view), so the two can never
  disagree and there is no copy-pasted logic.
- **`schema.sql` holds the DDL**, not Python string literals — the schema is easy
  to read and review on its own.
- **The SQL lives next to the question it answers** (window function for Part 2,
  cohort filter for Part 3, baseline join for Part 4), so each module reads like
  the question it implements.

---

## Database schema

Each CSV row is one sample. We normalise that flat row into four tables:

| Table | Grain | Key columns |
|---|---|---|
| `projects` | one project | `project_id` |
| `subjects` | one trial participant | `subject_id`, `project_id`, `condition`, `age`, `sex`, `treatment`, `response` |
| `samples` | one sample at one timepoint | `sample_id`, `subject_id`, `sample_type`, `time_from_treatment_start` |
| `cell_counts` | one (sample, population) | `sample_id`, `population`, `count` |

```
projects (1) ──< subjects (1) ──< samples (1) ──< cell_counts
```

### Design rationale

- **Third normal form.** Demographics, treatment, and response are *subject-level*
  facts — the EDA confirmed they never vary across a subject's samples — so they
  are stored **once per subject** instead of being repeated on every sample row.
  Updating a subject's response is a single-row write.
- **Counts are stored long, not wide.** `cell_counts` has one row per
  `(sample, population)` rather than five fixed columns. Adding a sixth population
  (or fifty) is *data*, not a schema migration — and computing per-sample totals or
  per-population statistics is a simple `GROUP BY`.
- **Foreign keys + a `CHECK (count >= 0)`** keep the data trustworthy; the loader
  inserts parents before children and `PRAGMA foreign_keys = ON` is set on every
  connection.
- **Indexes** on the columns the cohort queries filter and group by
  (`subjects.project_id`; `subjects(condition, treatment, response)`;
  `samples.subject_id`; `samples(sample_type, time_from_treatment_start)`;
  `cell_counts.population`).

### How this scales

With **hundreds of projects / thousands of samples**, this layout holds up:

- The normalised + indexed design keeps the Part 3/4 filter-and-group queries fast
  as `cell_counts` grows into the millions of rows; the database engine uses the
  indexes instead of scanning.
- The long `cell_counts` table absorbs **new assays/populations as rows**, so new
  analytics don't require schema changes or wide, sparse tables.
- Clean foreign keys let analysts join with confidence and let the data load
  incrementally per project/batch.
- If volume eventually outgrows SQLite, the identical schema and SQL port directly
  to Postgres/DuckDB; only the connection in `analysis/db.py` changes.

---

## Results summary

**Part 2** — `outputs/cell_frequencies.csv`: columns `sample, total_count,
population, count, percentage` (one row per sample × population). Each sample's
five percentages sum to 100.

**Part 3** — cohort = melanoma **and** miraclib **and** PBMC = **1,968 samples**
(993 responders, 975 non-responders). Mann-Whitney U test per population (a
non-parametric rank test — cell fractions are bounded and not guaranteed normal):

| population | median R | median NR | p (raw) | p (Bonferroni) |
|---|---|---|---|---|
| b_cell | 9.43 | 9.79 | 0.056 | 0.279 |
| cd8_t_cell | 24.73 | 24.60 | 0.639 | 1.000 |
| **cd4_t_cell** | **30.22** | **29.66** | **0.013** | 0.067 |
| nk_cell | 14.51 | 14.80 | 0.121 | 0.605 |
| monocyte | 19.61 | 19.94 | 0.163 | 0.816 |

**`cd4_t_cell` is the population that differs** between responders and
non-responders, with responders trending higher. Two honest caveats, both reported
in `outputs/`:

1. **Multiple testing.** We run five tests; under a conservative Bonferroni
   correction (α = 0.01) cd4_t_cell's adjusted p = 0.067 is **not** significant.
   The result is a *nominal* signal worth following up, not a confirmed biomarker.
2. **Repeated measures.** Each subject contributes several timepoints, so the
   sample-level test treats non-independent samples as independent. The
   subject-level sensitivity analysis (`part3_stats_subject_level.csv`) averages
   each subject's samples first; cd4_t_cell stays nominally significant
   (p ≈ 0.012) but again does not survive Bonferroni.

See `outputs/part3_boxplots.png` for the responder-vs-non-responder boxplots.

**Part 4** — baseline subset (melanoma, miraclib, PBMC, `time_from_treatment_start = 0`)
= **656 samples**:

- Samples per project: **prj1 = 384, prj3 = 272** (`part4_samples_per_project.csv`).
  prj2 contributes none — it has no melanoma + miraclib PBMC samples at baseline,
  so only prj1 and prj3 appear.
- Subjects by response: **yes = 331, no = 325** (`part4_subjects_by_response.csv`)
- Subjects by sex: **M = 344, F = 312** (`part4_subjects_by_sex.csv`)
- **Average B cells for melanoma male responders at baseline = 10401.28**
  (`part4_male_responder_bcells.csv`)

---

## Exploratory analysis

`datanalysis/` contains the exploration that drove the schema design
(`explore.py` → `eda_report.txt`, `fig_part3_boxplots.png`) and a script that
verifies the statistical claims (`verify_codex_stats.py`). These are documentation
of the design process; the authoritative outputs come from `make pipeline`.
