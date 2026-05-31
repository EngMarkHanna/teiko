"""Verify the statistical claims added in the Codex review (Section 9)."""
import os
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "cell-count.csv"))
POPS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
NTESTS = len(POPS)

long = df.melt(id_vars=["sample", "subject", "condition", "treatment",
                        "sample_type", "response", "time_from_treatment_start"],
               value_vars=POPS, var_name="population", value_name="count")
long["total"] = long.groupby("sample")["count"].transform("sum")
long["pct"] = long["count"] / long["total"] * 100

cohort = long[(long.condition == "melanoma") & (long.treatment == "miraclib")
              & (long.sample_type == "PBMC") & (long.response.isin(["yes", "no"]))]

def run(frame, label):
    print(f"\n=== {label} ===")
    print(f"{'population':>12} {'p_raw':>10} {'p_bonf':>10} {'sig0.05':>8} {'sigBonf':>8}")
    for pop in POPS:
        R = frame[(frame.population == pop) & (frame.response == "yes")].pct
        NR = frame[(frame.population == pop) & (frame.response == "no")].pct
        _, p = stats.mannwhitneyu(R, NR, alternative="two-sided")
        padj = min(p * NTESTS, 1.0)
        print(f"{pop:>12} {p:10.4f} {padj:10.4f} "
              f"{'YES' if p < 0.05 else 'no':>8} {'YES' if padj < 0.05 else 'no':>8}")

# (1) sample-level (what we already reported)
run(cohort, "sample-level (n samples R/NR = "
    f"{cohort[cohort.response=='yes']['sample'].nunique()}/"
    f"{cohort[cohort.response=='no']['sample'].nunique()})")

# (2) subject-level sensitivity: average each subject's pct across PBMC timepoints
subj = (cohort.groupby(["subject", "response", "population"], as_index=False)
              .pct.mean())
print(f"\nsubject-level cohort size: {subj.subject.nunique()} subjects "
      f"(R={subj[(subj.response=='yes')&(subj.population=='b_cell')].shape[0]}, "
      f"NR={subj[(subj.response=='no')&(subj.population=='b_cell')].shape[0]})")
run(subj, "subject-level (mean per subject)")

# (3) baseline-only PBMC
base = cohort[cohort.time_from_treatment_start == 0]
run(base, "baseline-only PBMC (time=0, n samples = "
    f"{base['sample'].nunique()})")
