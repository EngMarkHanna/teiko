"""
run_pipeline.py  --  run Parts 2-4 end to end and write outputs/.

    python run_pipeline.py

If the database is missing it is built first (Part 1), so this script alone
reproduces every table and plot. It finishes with regression checks that fail
loudly if any verified number drifts.
"""
import sys

import pandas as pd

import load_data
from analysis import frequencies, stats_tests, subsets
from analysis.db import DB_PATH, OUTPUTS_DIR


def _ensure_database() -> None:
    if not DB_PATH.exists():
        print("Database not found -> running load_data first.\n")
        if load_data.main() != 0:
            raise SystemExit("load_data failed")
        print()


def part2() -> pd.DataFrame:
    print("== Part 2: relative frequency table ==")
    freq = frequencies.cell_frequencies()
    out = OUTPUTS_DIR / "cell_frequencies.csv"
    freq.to_csv(out, index=False)
    print(f"  {len(freq):,} rows -> {out.name}")
    print(freq.head(5).to_string(index=False))
    return freq


def part3():
    print("\n== Part 3: responders vs non-responders ==")
    cohort = stats_tests.responder_cohort()
    sample_stats = stats_tests.compare_responders(cohort)
    subject_stats = stats_tests.compare_responders_by_subject(cohort)

    sample_stats.to_csv(OUTPUTS_DIR / "part3_stats_sample_level.csv", index=False)
    subject_stats.to_csv(OUTPUTS_DIR / "part3_stats_subject_level.csv", index=False)
    stats_tests.save_boxplot(cohort, OUTPUTS_DIR / "part3_boxplots.png")

    print(f"  cohort: {cohort['sample'].nunique():,} samples")
    print("  sample-level Mann-Whitney U:")
    print(sample_stats[["population", "median_responder", "median_non_responder",
                        "p_value", "p_bonferroni", "significant_p05",
                        "significant_bonferroni"]].to_string(index=False))
    raw_sig = stats_tests.significant_populations(sample_stats)
    bonf_sig = stats_tests.significant_populations(sample_stats, bonferroni=True)
    print(f"  significant (raw p<0.05):   {raw_sig or 'none'}")
    print(f"  significant (Bonferroni):   {bonf_sig or 'none'}")
    print("  saved: part3_boxplots.png, part3_stats_sample_level.csv, "
          "part3_stats_subject_level.csv")
    return cohort, sample_stats


def part4():
    print("\n== Part 4: baseline subset ==")
    subset = subsets.baseline_subset()
    per_project = subsets.samples_per_project(subset)
    by_response = subsets.subjects_by_response(subset)
    by_sex = subsets.subjects_by_sex(subset)
    male_resp_bcells = subsets.male_responder_mean_bcells(subset)

    per_project.to_csv(OUTPUTS_DIR / "part4_samples_per_project.csv", index=False)
    by_response.to_csv(OUTPUTS_DIR / "part4_subjects_by_response.csv", index=False)
    by_sex.to_csv(OUTPUTS_DIR / "part4_subjects_by_sex.csv", index=False)
    pd.DataFrame([{
        "metric": "avg_b_cells_melanoma_male_responders_baseline",
        "value": male_resp_bcells,
    }]).to_csv(OUTPUTS_DIR / "part4_male_responder_bcells.csv", index=False)

    print(f"  subset: {len(subset):,} samples")
    print("  samples per project:\n" + per_project.to_string(index=False))
    print("  subjects by response:\n" + by_response.to_string(index=False))
    print("  subjects by sex:\n" + by_sex.to_string(index=False))
    print(f"  avg B cells, melanoma male responders at baseline: {male_resp_bcells:.2f}")
    return subset


def regression_checks(freq, cohort, sample_stats, subset) -> None:
    """Fail loudly if any verified, deterministic number has drifted."""
    print("\n== Regression checks ==")

    # Part 2
    assert len(freq) == 52_500, f"freq rows {len(freq)} != 52500"
    pct_sums = freq.groupby("sample").percentage.sum().round(2)
    assert (pct_sums == 100.00).all(), "per-sample percentages do not sum to 100"

    # Part 3
    by_resp = cohort.groupby("response")["sample"].nunique()
    assert cohort["sample"].nunique() == 1968, "cohort != 1968 samples"
    assert by_resp.get("yes") == 993, "responders != 993"
    assert by_resp.get("no") == 975, "non-responders != 975"
    assert stats_tests.significant_populations(sample_stats) == ["cd4_t_cell"], \
        "raw-significant populations changed"

    # Part 4
    assert len(subset) == 656, f"baseline subset {len(subset)} != 656"
    proj = subset.groupby("project_id").size()
    assert proj.get("prj1") == 384, "prj1 != 384"
    assert proj.get("prj3") == 272, "prj3 != 272"
    avg = subsets.male_responder_mean_bcells(subset)
    assert abs(avg - 10401.28) < 0.005, f"male-responder B-cell avg {avg} != 10401.28"

    print("  all checks passed.")


def main() -> int:
    _ensure_database()
    OUTPUTS_DIR.mkdir(exist_ok=True)

    freq = part2()
    cohort, sample_stats = part3()
    subset = part4()
    regression_checks(freq, cohort, sample_stats, subset)

    print(f"\nDone. Outputs written to {OUTPUTS_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
