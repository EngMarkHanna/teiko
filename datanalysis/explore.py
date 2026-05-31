"""
Exploratory data analysis for cell-count.csv.

Goal: understand the dataset *before* designing the schema / pipeline, and
pre-compute the deterministic answers (Parts 2-4) so we can verify the real
pipeline later. Uses pandas / scipy / seaborn from the project venv.

Run:  python datanalysis/explore.py
Outputs:
    datanalysis/eda_report.txt           human-readable findings
    datanalysis/fig_part3_boxplots.png   responder vs non-responder boxplots
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "cell-count.csv")
REPORT = os.path.join(HERE, "eda_report.txt")
FIG = os.path.join(HERE, "fig_part3_boxplots.png")

POPS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

_log = []
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    _log.append(s)

df = pd.read_csv(CSV)

# ---------------------------------------------------------------- shape
out("=" * 72)
out("SHAPE & COLUMNS")
out("=" * 72)
out("rows:", len(df), "| cols:", df.shape[1])
out("columns:", list(df.columns))
out("\ndtypes:\n" + df.dtypes.to_string())

# ---------------------------------------------------------- categoricals
out("\n" + "=" * 72)
out("CATEGORICAL VALUE COUNTS")
out("=" * 72)
for c in ["project", "condition", "sex", "treatment", "response",
          "sample_type", "time_from_treatment_start"]:
    out(f"\n[{c}]")
    out(df[c].value_counts(dropna=False).to_string())

# --------------------------------------------------------------- missing
out("\n" + "=" * 72)
out("MISSING VALUES (only non-zero shown)")
out("=" * 72)
miss = df.isna().sum()
out(miss[miss > 0].to_string() if (miss > 0).any() else "none")
out("\nNote: 'response' missing is expected for non-treated / healthy subjects.")

# ----------------------------------------------------- relational checks
out("\n" + "=" * 72)
out("RELATIONAL INTEGRITY (supports the 4-table schema)")
out("=" * 72)
out("unique projects:", df.project.nunique())
out("unique subjects:", df.subject.nunique())
out("unique samples :", df["sample"].nunique(), "| total rows:", len(df),
    "->", "sample is unique" if df["sample"].is_unique else "DUPLICATE samples!")

subj_proj = df.groupby("subject").project.nunique()
out("subjects spanning >1 project:", int((subj_proj > 1).sum()))
attr_cols = ["condition", "age", "sex", "treatment", "response"]
inconsistent = df.groupby("subject")[attr_cols].nunique().gt(1).any(axis=1).sum()
out("subjects whose demographic attrs vary across their samples:", int(inconsistent))
out("=> subject-level columns can live in a `subjects` table without duplication"
    if inconsistent == 0 else "=> attrs vary; revisit schema")

neg = (df[POPS] < 0).any().any()
out("any negative cell counts:", bool(neg))

# ----------------------------------------------------- PART 2 long table
out("\n" + "=" * 72)
out("PART 2  -  relative frequency summary (sample/total_count/population/count/percentage)")
out("=" * 72)
long = df.melt(id_vars=["sample"], value_vars=POPS,
               var_name="population", value_name="count")
long["total_count"] = long.groupby("sample")["count"].transform("sum")
long["percentage"] = (long["count"] / long["total_count"] * 100).round(4)
long = long[["sample", "total_count", "population", "count", "percentage"]]
out("rows in long table:", len(long), "(=", len(df), "samples x 5 pops)")
out("\nhead:\n" + long.head(10).to_string(index=False))
sums = long.groupby("sample").percentage.sum().round(2)
out("\nper-sample percentage sum  -> min:", sums.min(), "max:", sums.max(),
    "(should be 100)")

# --------------------------------------------------- PART 3 statistics
out("\n" + "=" * 72)
out("PART 3  -  melanoma & miraclib & PBMC : responders vs non-responders")
out("=" * 72)
meta = df[["sample", "condition", "treatment", "sample_type", "response",
           "project", "subject", "sex", "time_from_treatment_start"]]
p3 = long.merge(meta, on="sample")
p3 = p3[(p3.condition == "melanoma") & (p3.treatment == "miraclib")
        & (p3.sample_type == "PBMC") & (p3.response.isin(["yes", "no"]))]
out("cohort samples:", p3["sample"].nunique())
out("  responder samples:", p3[p3.response == "yes"]["sample"].nunique())
out("  non-responder samples:", p3[p3.response == "no"]["sample"].nunique())

out("\nMann-Whitney U (two-sided) per population, on percentage:")
out(f"  {'population':>12} {'med_R':>8} {'med_NR':>8} {'U':>10} {'p_value':>12} {'sig?':>6}")
results = []
for pop in POPS:
    R = p3[(p3.population == pop) & (p3.response == "yes")].percentage
    NR = p3[(p3.population == pop) & (p3.response == "no")].percentage
    U, pval = stats.mannwhitneyu(R, NR, alternative="two-sided")
    sig = "YES" if pval < 0.05 else "no"
    results.append((pop, R.median(), NR.median(), U, pval, sig))
    out(f"  {pop:>12} {R.median():8.2f} {NR.median():8.2f} {U:10.1f} {pval:12.3e} {sig:>6}")
sig_pops = [r[0] for r in results if r[5] == "YES"]
out("\nSIGNIFICANT (p<0.05):", sig_pops if sig_pops else "none")

# boxplots
plt.figure(figsize=(12, 6))
sns.boxplot(data=p3, x="population", y="percentage", hue="response",
            order=POPS, hue_order=["yes", "no"])
plt.title("Relative frequency by population: responders (yes) vs non-responders (no)\n"
          "melanoma / miraclib / PBMC")
plt.ylabel("relative frequency (%)")
plt.xlabel("immune cell population")
plt.tight_layout()
plt.savefig(FIG, dpi=120)
plt.close()
out("\n[boxplot saved -> datanalysis/fig_part3_boxplots.png]")

# ----------------------------------------------------- PART 4 subset
out("\n" + "=" * 72)
out("PART 4  -  melanoma PBMC baseline (time=0) miraclib")
out("=" * 72)
p4 = df[(df.condition == "melanoma") & (df.treatment == "miraclib")
        & (df.sample_type == "PBMC") & (df.time_from_treatment_start == 0)]
out("samples in subset:", len(p4))

out("\nsamples per project:")
out(p4.project.value_counts().sort_index().to_string())

subj = p4.drop_duplicates("subject")
out("\nunique subjects:", len(subj))
out("subjects by response:")
out(subj.response.value_counts(dropna=False).to_string())
out("subjects by sex:")
out(subj.sex.value_counts(dropna=False).to_string())

mr = p4[(p4.sex == "M") & (p4.response == "yes")]
out("\nmelanoma MALES, responders, t=0  -> n samples:", len(mr))
out("AVERAGE B CELLS (XXX.XX): {:.2f}".format(mr.b_cell.mean()))

# ----------------------------------------------------- write report
with open(REPORT, "w") as f:
    f.write("\n".join(_log) + "\n")
out("\n[report saved -> datanalysis/eda_report.txt]")
