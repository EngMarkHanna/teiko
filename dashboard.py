"""
dashboard.py  --  interactive Streamlit dashboard for Parts 2-4.

    streamlit run dashboard.py     (or: make dashboard)

Reuses the exact same analysis functions as run_pipeline.py, so the dashboard
and the batch outputs can never disagree. If the database is missing it is
built on first load.
"""
import matplotlib.pyplot as plt
import streamlit as st

import load_data
from analysis import frequencies, stats_tests, subsets
from analysis.db import DB_PATH, POPULATIONS

st.set_page_config(page_title="Loblaw Bio - Cell Count Dashboard", layout="wide")


# Cached loaders so the widgets stay snappy.
@st.cache_resource
def ensure_database():
    if not DB_PATH.exists():
        load_data.main()
    return True


@st.cache_data
def load_frequencies():
    return frequencies.cell_frequencies()


@st.cache_data
def load_cohort():
    return stats_tests.responder_cohort()


@st.cache_data
def load_subset():
    return subsets.baseline_subset()


ensure_database()

st.title("Loblaw Bio - Immune Cell Population Dashboard")
st.caption("Relative cell-population frequencies, responder analysis, and the "
           "baseline subset, computed live from the SQLite database.")

tab2, tab3, tab4 = st.tabs([
    "Part 2 - Frequencies",
    "Part 3 - Responders vs Non-responders",
    "Part 4 - Baseline subset",
])

with tab2:
    st.subheader("Relative frequency of each population per sample")
    freq = load_frequencies()

    c1, c2 = st.columns(2)
    c1.metric("Samples", f"{freq['sample'].nunique():,}")
    c2.metric("Rows (sample x population)", f"{len(freq):,}")

    pops = st.multiselect("Populations", POPULATIONS, default=POPULATIONS)
    search = st.text_input("Filter by sample id (substring)", "")
    view = freq[freq.population.isin(pops)]
    if search:
        view = view[view["sample"].str.contains(search, case=False)]
    st.dataframe(view, width="stretch", hide_index=True)

    st.markdown("**Composition of a single sample**")
    sample_id = st.selectbox("Sample", sorted(freq["sample"].unique()))
    one = (freq[freq["sample"] == sample_id]
           .set_index("population")["percentage"]
           .reindex(POPULATIONS))
    st.bar_chart(one)

    st.markdown("**Compare samples per population (grouped bars)**")
    all_samples = sorted(freq["sample"].unique())
    max_compare = 10
    chosen = st.multiselect(
        f"Pick up to {max_compare} samples to compare side by side",
        options=all_samples,
        default=all_samples[:5],
        max_selections=max_compare,
    )
    if chosen:
        # Group by population on the x-axis, one bar per sample within each group.
        grouped = (freq[freq["sample"].isin(chosen)]
                   .pivot(index="population", columns="sample", values="percentage")
                   .reindex(index=POPULATIONS, columns=chosen))
        st.bar_chart(grouped, stack=False,
                     x_label="immune cell population",
                     y_label="relative frequency (%)")
        st.caption("Each population is a group on the x-axis; within it, one bar "
                   "per selected sample, so you can compare samples zone by zone.")
    else:
        st.info("Select at least one sample to see the comparison.")

with tab3:
    st.subheader("Melanoma / miraclib / PBMC: responders vs non-responders")
    cohort = load_cohort()
    sample_stats = stats_tests.compare_responders(cohort)
    subject_stats = stats_tests.compare_responders_by_subject(cohort)

    by_resp = cohort.groupby("response")["sample"].nunique()
    c1, c2, c3 = st.columns(3)
    c1.metric("Cohort samples", f"{cohort['sample'].nunique():,}")
    c2.metric("Responder samples", f"{int(by_resp.get('yes', 0)):,}")
    c3.metric("Non-responder samples", f"{int(by_resp.get('no', 0)):,}")

    fig, ax = plt.subplots(figsize=(11, 5))
    stats_tests.boxplot(cohort, ax=ax)
    st.pyplot(fig)

    st.markdown("**Mann-Whitney U test (sample level)**")
    st.dataframe(sample_stats, width="stretch", hide_index=True)

    raw_sig = stats_tests.significant_populations(sample_stats)
    bonf_sig = stats_tests.significant_populations(sample_stats, bonferroni=True)
    if raw_sig:
        st.success(f"Significant at raw p<0.05: **{', '.join(raw_sig)}**")
    if not bonf_sig:
        st.warning("No population survives Bonferroni correction (5 tests, "
                   "alpha=0.01) - treat the result as nominal, not confirmatory.")

    with st.expander("Subject-level sensitivity (average each subject's timepoints)"):
        st.dataframe(subject_stats, width="stretch", hide_index=True)
        st.caption("Samples are repeated measures per subject. Averaging each "
                   "subject first gives independent observations; cd4_t_cell stays "
                   "nominally significant but still does not survive Bonferroni.")

with tab4:
    st.subheader("Melanoma PBMC baseline (time = 0), miraclib-treated")
    subset = load_subset()

    c1, c2 = st.columns(2)
    c1.metric("Samples in subset", f"{len(subset):,}")
    c2.metric("Avg B cells - male responders",
              f"{subsets.male_responder_mean_bcells(subset):.2f}")

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**Samples per project**")
        st.bar_chart(subsets.samples_per_project(subset)
                     .set_index("project_id")["n_samples"])
    with c2:
        st.markdown("**Subjects by response**")
        st.bar_chart(subsets.subjects_by_response(subset)
                     .set_index("response")["n_subjects"])
    with c3:
        st.markdown("**Subjects by sex**")
        st.bar_chart(subsets.subjects_by_sex(subset)
                     .set_index("sex")["n_subjects"])

    st.caption("Note: prj2 contributes no samples to this subset — it has no "
               "melanoma + miraclib PBMC samples at baseline (time = 0), so only "
               "prj1 and prj3 appear above.")
