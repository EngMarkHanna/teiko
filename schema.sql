-- Relational model for the cell-count dataset: one flat CSV row becomes four
-- tables so subject-level facts are stored once and counts stay long/tidy.
-- Applied by load_data.py.

-- Drop child -> parent order so foreign keys never dangle.
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id  TEXT PRIMARY KEY
);

-- Demographics/treatment/response are subject-level (constant across a
-- subject's samples), so they live here once rather than on every sample.
CREATE TABLE subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL REFERENCES projects(project_id),
    condition   TEXT,            -- melanoma / carcinoma / healthy
    age         INTEGER,
    sex         TEXT,            -- M / F
    treatment   TEXT,            -- miraclib / phauximab / none
    response    TEXT             -- 'yes' / 'no' / NULL (untreated/healthy)
);

-- One sample drawn from a subject at one timepoint.
CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT,        -- PBMC / WB
    time_from_treatment_start  INTEGER      -- 0 / 7 / 14
);

-- Long/tidy: one row per (sample, population), so a new population is data,
-- not a schema change.
CREATE TABLE cell_counts (
    sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
    population  TEXT NOT NULL,   -- b_cell / cd8_t_cell / cd4_t_cell / nk_cell / monocyte
    count       INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

-- Indexes on the columns the Part 3/4 cohort queries filter and group by.
CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_condition_treatment_response
    ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_samples_type_time ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_population ON cell_counts(population);
