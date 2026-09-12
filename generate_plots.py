#!/usr/bin/env python3
"""Generate figures from the CSVs under work/*/results.

The script has two modes, because the figures serve two different purposes.

``--mode report`` (the default) builds the figures for the write-up: it reads
only the full-size ``paper_subset_50x10_*`` designs, which partition cleanly by
sequence length (one design per length), so rows are tagged with their sequence
length and every figure facets on that rather than on the design name. It plots
per-replicate rows with 95% CIs, uses typeset model names, and writes vector
PDFs to ``report/figures``.

``--mode explore`` is the survey view: it reads every design directory that has
a ``results/`` folder, tags rows with the design name, facets on design, and
writes PNGs to ``work/plots`` -- plus the cross-design summaries (sequence-length
trend, model ranking, pairwise error matrix) that only make sense once several
designs are in play.

Usage:
    .venv-plots/bin/python generate_plots.py
    .venv-plots/bin/python generate_plots.py --mode explore --split-by-design

--- Tweakable parameters ---
Everything that controls *which* models/designs/ILS levels/sequence lengths go
into the plots lives in the CONFIG block right below the imports. Edit those
lists directly, or override them for a single run with the matching ``--models``
/ ``--designs`` / ``--exclude-designs`` / ``--ils`` / ``--seq-lengths`` CLI flags
(comma-separated). Leaving a CLI flag unset keeps the CONFIG default. The
filters apply in both modes. ``--split-by-design`` (or SPLIT_BY_DESIGN) is
explore-only: it writes one PNG per design instead of combining all designs
into one figure.
"""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ============================= CONFIG =====================================
# All models known to the pipeline, in the order they should appear on axes
# and legends. Trim this list to drop a model from every plot.
MODEL_ORDER = ["gtr_g4", "gtr", "hky_g4", "hky", "jc"]

# Which of the models above are actually included in a given run. Defaults
# to "all of MODEL_ORDER"; set to e.g. ["gtr_g4", "gtr"] to only compare
# those two everywhere.
INCLUDE_MODELS = list(MODEL_ORDER)

# Which design directories under --work-dir to include/exclude. Matched as
# substrings against the directory name (e.g. "paper_subset_50x10" matches
# "paper_subset_50x10_L200", "paper_subset_50x10_L800", ...).
# None means "include everything found on disk".
INCLUDE_DESIGNS: list[str] | None = None
EXCLUDE_DESIGNS: list[str] = []

# Report mode overrides INCLUDE_DESIGNS with this: the full-size runs are the
# only ones that carry one sequence length per design, which is the layout
# every report figure facets on.
REPORT_DESIGNS = ["paper_subset_50x10"]

# Which ILS levels ("low", "high") to keep. None means "keep all".
INCLUDE_ILS: list[str] | None = None

# Which sequence lengths (as they appear in the "sequence_length" column) to
# keep. None means "keep all". Only applies to CSVs that have that column.
INCLUDE_SEQ_LENGTHS: list[int] | None = None

# By default, plots that facet by design (e.g. stage1_high_support_error,
# stage2_species_tree_error) cram every included design into one PNG as
# subplot rows/columns. Set this to True (or pass --split-by-design) to
# instead write one separate PNG per design, suffixed "__<design name>".
# Explore mode only.
SPLIT_BY_DESIGN = False

# Typeset names for the report figures; explore mode keeps the raw column
# values, which are what the CSV filters and filenames use.
MODEL_LABELS = {
    "gtr_g4": "GTR+Γ4",
    "gtr": "GTR",
    "hky_g4": "HKY+Γ4",
    "hky": "HKY",
    "jc": "JC",
}
ILS_ORDER = ["low", "high"]

DEFAULT_OUT = {"report": Path("report/figures"), "explore": Path("work/plots")}
# ===========================================================================

sns.set_theme(style="whitegrid", context="talk")
matplotlib.rcParams["pdf.fonttype"] = 42


# =========================== shared loading ================================


def _matches_any(name: str, substrings: list[str]) -> bool:
    return any(s in name for s in substrings)


def load_all(work_dir: Path, csv_name: str,
             include_designs: list[str] | None = None,
             exclude_designs: list[str] | None = None) -> pd.DataFrame:
    frames = []
    for design_dir in sorted(work_dir.iterdir()):
        if not design_dir.is_dir():
            continue
        name = design_dir.name
        if include_designs and not _matches_any(name, include_designs):
            continue
        if exclude_designs and _matches_any(name, exclude_designs):
            continue
        csv_path = design_dir / "results" / csv_name
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        df.insert(0, "design", name)
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No {csv_name} found under {work_dir}/*/results")
    return pd.concat(frames, ignore_index=True)


def apply_filters(df: pd.DataFrame, models: list[str] | None = None,
                   ils: list[str] | None = None,
                   seq_lengths: list[int] | None = None) -> pd.DataFrame:
    """Filter rows to the configured subset of models/ILS/sequence lengths."""
    df = df.copy()
    if models is not None and "analysis_model" in df.columns:
        df = df[df["analysis_model"].isin(models)]
    if ils is not None and "ils" in df.columns:
        df = df[df["ils"].isin(ils)]
    if seq_lengths is not None and "sequence_length" in df.columns:
        df = df[df["sequence_length"].isin(seq_lengths)]
    return df


def order_models(df: pd.DataFrame, col: str = "analysis_model",
                  model_order: list[str] | None = None) -> pd.DataFrame:
    order = model_order or MODEL_ORDER
    df[col] = pd.Categorical(df[col], categories=order, ordered=True)
    return df.sort_values(col)


def design_group(design_name: str) -> str:
    """Strip the trailing "_L<seq lengths>" suffix, e.g.

    "paper_subset_50x10_L1600" -> "paper_subset_50x10". Useful for comparing
    the same taxa/loci design across its different sequence lengths.
    """
    return re.sub(r"_L[\d_]+$", "", design_name)


def split_frames(dfs: list[pd.DataFrame], split: bool,
                  group_col: str = "design") -> list[tuple[str, list[pd.DataFrame]]]:
    """Yield (filename_suffix, filtered_dfs) pairs.

    When ``split`` is False, yields a single ("", dfs) pair unchanged (today's
    behaviour: one combined multi-panel figure per plot). When True, yields
    one pair per distinct value of ``group_col`` found across ``dfs``, each
    filtered to just that value, so callers can write one file per design
    instead of cramming every design into one figure.
    """
    if not split:
        return [("", dfs)]
    values: set[str] = set()
    for df in dfs:
        if group_col in df.columns:
            values.update(df[group_col].unique())
    out = []
    for value in sorted(values):
        filtered = [df[df[group_col] == value].copy() if group_col in df.columns else df.copy()
                    for df in dfs]
        out.append((f"__{value}", filtered))
    return out


def parse_calibration_condition(cal: pd.DataFrame) -> pd.DataFrame:
    """Split the packed ``condition`` string into columns of its own.

    Report mode needs these before filtering, since ``stage1_calibration.csv``
    encodes model/ILS/length into one string rather than as real columns.
    """
    df = cal.copy()
    df["analysis_model"] = df["condition"].str.extract(r"model=([a-z_0-9]+)")
    df["ils"] = df["condition"].str.extract(r"ils=([a-z]+)")
    df["sequence_length"] = df["condition"].str.extract(r"length=(\d+)").astype(int)
    return df


# ===================== explore mode: PNGs, faceted by design ===============


def safe_aspect(n_cols: int, height: float, base_aspect: float, min_width: float = 6.0) -> float:
    """Widen narrow FacetGrids (e.g. only one ``col`` value left after filtering)

    so seaborn's "legend outside the grid" placement always has enough room;
    with a too-narrow figure it can raise "left cannot be >= right" from
    matplotlib when computing the legend's reserved margin.
    """
    n_cols = max(n_cols, 1)
    needed = min_width / (n_cols * height)
    return max(base_aspect, needed)


def finalize_facetgrid(g: sns.FacetGrid, out_dir: Path, name: str,
                        ylabel: str, xlabel: str,
                        legend_right: bool = False) -> None:
    """Apply consistent, collision-free layout to a row+col FacetGrid and save.

    Fixes three problems the original plots all shared:
    - long per-axes y-axis labels repeated on every row, overlapping each
      other vertically -> replaced with a single figure-level supylabel.
    - facet titles colliding between adjacent columns -> shortened via
      set_titles (done by the caller) plus rotated x tick labels so model
      names don't run into each other.
    - "outside" legends drawn on top of an interior panel with catplot's
      legend_out under row+col faceting -> the figure is widened and the
      legend is explicitly re-anchored in figure-fraction coordinates.
    """
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelrotation=30)
        ax.title.set_size(11)
        ax.set_ylabel("")
        ax.set_xlabel("")
    g.figure.supylabel(ylabel)
    g.figure.supxlabel(xlabel)

    fig_w, fig_h = g.figure.get_size_inches()
    right_margin = 0.8
    if legend_right and g._legend is not None:
        # Long legend labels ("support-only wASTRAL - contracted ASTRAL")
        # need real room. Rather than guessing a width from character count
        # (which undershot and still clipped the text), measure the legend's
        # actual rendered extent after a draw pass and size the figure to it.
        g.figure.canvas.draw()
        legend_bbox = g._legend.get_window_extent()
        legend_w_in = legend_bbox.width / g.figure.dpi
        margin_in = 0.35
        g.figure.set_size_inches(fig_w + legend_w_in + margin_in, fig_h + 0.5)
        right_margin = fig_w / (fig_w + legend_w_in + margin_in)
        g._legend.set_bbox_to_anchor((right_margin, 0.6), transform=g.figure.transFigure)
        g._legend.set_loc("center left")
    else:
        g.figure.set_size_inches(fig_w, fig_h + 0.5)
        right_margin = 0.98
    g.figure.subplots_adjust(top=0.90, right=right_margin, left=0.1, bottom=0.12,
                              hspace=0.55, wspace=0.3)

    path = out_dir / f"{name}.png"
    g.figure.savefig(path, dpi=150)
    plt.close(g.figure)
    print(f"wrote {path}")


def save(fig, out_dir: Path, name: str) -> None:
    path = out_dir / f"{name}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


def plot_branch_error(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                       suffix: str = "") -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy(), model_order=model_order)
    g = sns.catplot(
        data=pooled, x="analysis_model", y="branch_error_rate", hue="ils",
        col="design", kind="bar", height=4.5, aspect=0.9, col_wrap=2,
        order=model_order, legend=True,
    )
    g.set_axis_labels("Analysis model", "Branch error rate")
    g.set_titles("{col_name}", size=11)
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelrotation=30)
    # Re-anchor the auto-generated legend outside the axes grid: with
    # col_wrap, seaborn's legend_out=True placed it on top of the second
    # panel rather than beside the whole figure, which is what caused the
    # collision in the original plot.
    fig_w, fig_h = g.figure.get_size_inches()
    g.figure.set_size_inches(fig_w + 1.4, fig_h + 0.4)
    g.figure.subplots_adjust(top=0.92, right=0.8, bottom=0.12, hspace=0.55, wspace=0.3)
    if g._legend is not None:
        g._legend.set_bbox_to_anchor((0.82, 0.95), transform=g.figure.transFigure)
        g._legend.set_loc("upper left")
    path = out_dir / f"stage1_branch_error_by_model{suffix}.png"
    g.figure.savefig(path, dpi=150)
    plt.close(g.figure)
    print(f"wrote {path}")


def plot_calibration_metrics(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                              suffix: str = "") -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy(), model_order=model_order)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sns.barplot(data=pooled, x="analysis_model", y="brier_score", hue="ils",
                order=model_order, ax=axes[0])
    axes[0].set_title("Brier score (lower = better)")
    axes[0].set_xlabel("Analysis model")
    sns.barplot(data=pooled, x="analysis_model", y="expected_calibration_error", hue="ils",
                order=model_order, ax=axes[1])
    axes[1].set_title("Expected calibration error")
    axes[1].set_xlabel("Analysis model")
    save(fig, out_dir, f"stage1_calibration_scores{suffix}")


def plot_high_support_wrong(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                             suffix: str = "") -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy(), model_order=model_order)
    melted = pooled.melt(
        id_vars=["design", "ils", "analysis_model"],
        value_vars=["p_wrong_given_support_ge_0.95", "p_wrong_given_support_ge_0.99"],
        var_name="threshold", value_name="p_wrong",
    )
    aspect = safe_aspect(melted["ils"].nunique(), 3.5, 1.2)
    g = sns.catplot(
        data=melted, x="analysis_model", y="p_wrong", hue="threshold",
        col="ils", row="design", kind="bar", order=model_order, height=3.5, aspect=aspect,
        legend_out=True,
    )
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, f"stage1_high_support_error{suffix}",
                        ylabel="P(wrong | support >= t)", xlabel="Analysis model",
                        legend_right=True)


def plot_calibration_curve(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                            suffix: str = "") -> None:
    df = df.copy()
    df["model"] = df["condition"].str.extract(r"model=([a-z_0-9]+)")
    df["ils"] = df["condition"].str.extract(r"ils=([a-z]+)")
    df = df[df["model"].isin(model_order)]
    df["model"] = pd.Categorical(df["model"], categories=model_order, ordered=True)
    aspect = safe_aspect(df["ils"].nunique(), 3.5, 1.2)
    g = sns.relplot(
        data=df.sort_values("mean_support"), x="mean_support", y="p_correct",
        hue="model", hue_order=model_order, col="ils", row="design",
        kind="line", marker="o", height=3.5, aspect=aspect, facet_kws={"sharex": True, "sharey": True},
    )
    for ax in g.axes.flat:
        ax.plot([0, 1], [0, 1], ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, f"stage1_calibration_curve{suffix}",
                        ylabel="Observed P(correct)", xlabel="Mean support in bin",
                        legend_right=True)


def plot_stage2_nrf(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                     suffix: str = "") -> None:
    df = order_models(df.copy(), model_order=model_order)
    melted = df.melt(
        id_vars=["design", "ils", "replicate", "sequence_length", "analysis_model"],
        value_vars=["unweighted_nrf", "contracted_unweighted_nrf", "weighted_nrf", "hybrid_weighted_nrf"],
        var_name="method", value_name="nrf",
    )
    method_labels = {
        "unweighted_nrf": "ASTRAL (resolved)",
        "contracted_unweighted_nrf": "ASTRAL (0.90-contracted)",
        "weighted_nrf": "wASTRAL (support-only)",
        "hybrid_weighted_nrf": "wASTRAL (hybrid)",
    }
    melted["method"] = melted["method"].map(method_labels)
    aspect = safe_aspect(melted["ils"].nunique(), 3.5, 1.3)
    g = sns.catplot(
        data=melted, x="analysis_model", y="nrf", hue="method",
        col="ils", row="design", kind="box", order=model_order, height=3.5, aspect=aspect,
        legend_out=True,
    )
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, f"stage2_species_tree_error{suffix}",
                        ylabel="Normalized RF to true species tree", xlabel="Analysis model",
                        legend_right=True)


def plot_stage2_deltas(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                        suffix: str = "") -> None:
    df = order_models(df.copy(), model_order=model_order)
    melted = df.melt(
        id_vars=["design", "ils", "replicate", "sequence_length", "analysis_model"],
        value_vars=["delta_weighted_minus_contracted", "delta_hybrid_minus_contracted"],
        var_name="contrast", value_name="delta_nrf",
    )
    contrast_labels = {
        "delta_weighted_minus_contracted": "support-only wASTRAL - contracted ASTRAL",
        "delta_hybrid_minus_contracted": "hybrid wASTRAL - contracted ASTRAL",
    }
    melted["contrast"] = melted["contrast"].map(contrast_labels)
    aspect = safe_aspect(melted["ils"].nunique(), 3.5, 1.3)
    g = sns.catplot(
        data=melted, x="analysis_model", y="delta_nrf", hue="contrast",
        col="ils", row="design", kind="box", order=model_order, height=3.5, aspect=aspect,
        legend_out=True,
    )
    for ax in g.axes.flat:
        ax.axhline(0, ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, f"stage2_primary_contrasts{suffix}",
                        ylabel="Delta nRF (negative favors wASTRAL)", xlabel="Analysis model",
                        legend_right=True)


def plot_ils_separation(df: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sns.barplot(data=df, x="design", y="mean_true_gene_species_nrf", hue="ils", ax=ax)
    ax.set_ylabel("Mean true gene/species nRF")
    ax.set_xlabel("Design")
    ax.tick_params(axis="x", rotation=20)
    save(fig, out_dir, "ils_separation_check")


def plot_seq_length_effect(stage1: pd.DataFrame, stage2: pd.DataFrame,
                            out_dir: Path, model_order: list[str]) -> None:
    """How each model's error trends as sequence length grows, per taxa/loci design."""
    for label, df, metric, ylabel, name in (
        ("stage1", stage1, "branch_error_rate", "Branch error rate", "seq_length_effect_stage1"),
        ("stage2", stage2, "contracted_unweighted_nrf", "Species-tree nRF (contracted ASTRAL)",
         "seq_length_effect_stage2"),
    ):
        pooled = df
        if "replicate" in pooled.columns and pooled["replicate"].dtype == object:
            pooled = pooled[pooled["replicate"] == "pooled"]
        pooled = pooled.copy()
        pooled["design_group"] = pooled["design"].map(design_group)
        if pooled["design_group"].nunique() < 1 or pooled["sequence_length"].nunique() < 2:
            print(f"skipping {name}: not enough sequence-length variation")
            continue
        pooled = order_models(pooled, model_order=model_order)
        aspect = safe_aspect(pooled["ils"].nunique(), 3.5, 1.2)
        g = sns.relplot(
            data=pooled.sort_values("sequence_length"), x="sequence_length", y=metric,
            hue="analysis_model", hue_order=model_order, col="ils", row="design_group",
            kind="line", marker="o", height=3.5, aspect=aspect,
            facet_kws={"sharex": False, "sharey": True},
        )
        g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
        finalize_facetgrid(g, out_dir, name, ylabel=ylabel,
                            xlabel="Sequence length", legend_right=True)


def plot_model_ranking(stage1: pd.DataFrame, stage2: pd.DataFrame,
                        out_dir: Path, model_order: list[str]) -> None:
    """Average rank of each model (1 = best) across every design/ILS condition."""
    pooled1 = stage1[stage1["replicate"] == "pooled"] if stage1["replicate"].dtype == object else stage1
    rank1 = pooled1.copy()
    rank1["rank"] = rank1.groupby(["design", "ils"])["branch_error_rate"].rank()
    rank1["metric"] = "Stage 1 branch error rate"

    rank2 = stage2.copy()
    rank2["rank"] = rank2.groupby(["design", "ils", "replicate"])["contracted_unweighted_nrf"].rank()
    rank2["metric"] = "Stage 2 species-tree nRF"

    combined = pd.concat([
        rank1[["analysis_model", "rank", "metric"]],
        rank2[["analysis_model", "rank", "metric"]],
    ], ignore_index=True)
    combined = order_models(combined, model_order=model_order)

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.barplot(data=combined, x="analysis_model", y="rank", hue="metric",
                order=model_order, ax=ax, errorbar="ci")
    ax.set_ylabel("Mean rank (1 = best)")
    ax.set_xlabel("Analysis model")
    ax.tick_params(axis="x", rotation=30)
    save(fig, out_dir, "model_ranking_summary")


def plot_pairwise_model_heatmap(stage1: pd.DataFrame, out_dir: Path,
                                 model_order: list[str]) -> None:
    """Pairwise mean-branch-error-rate difference matrix between every model pair."""
    pooled = stage1[stage1["replicate"] == "pooled"] if stage1["replicate"].dtype == object else stage1
    means = pooled.groupby("analysis_model")["branch_error_rate"].mean().reindex(model_order)
    diff = pd.DataFrame(
        [[a - b for b in means] for a in means], index=model_order, columns=model_order,
    )
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(diff, annot=True, fmt=".3f", cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "Row minus column (branch error rate)"})
    ax.set_xlabel("Model (column)")
    ax.set_ylabel("Model (row)")
    save(fig, out_dir, "pairwise_model_branch_error_diff")


def plot_replicate_variability(df: pd.DataFrame, out_dir: Path, model_order: list[str],
                                suffix: str = "") -> None:
    """Spread of branch error rate across individual replicates (not pooled)."""
    per_rep = df[df["replicate"] != "pooled"] if df["replicate"].dtype == object else df
    if per_rep.empty:
        print("skipping stage1_replicate_variability: no per-replicate rows found")
        return
    per_rep = order_models(per_rep.copy(), model_order=model_order)
    aspect = safe_aspect(min(per_rep["design"].nunique(), 2), 4, 1.0)
    g = sns.catplot(
        data=per_rep, x="analysis_model", y="branch_error_rate", hue="ils",
        col="design", kind="box", order=model_order, height=4, aspect=aspect, col_wrap=2,
        legend_out=True,
    )
    g.set_titles("{col_name}", size=11)
    finalize_facetgrid(g, out_dir, f"stage1_replicate_variability{suffix}",
                        ylabel="Branch error rate (per replicate)", xlabel="Analysis model",
                        legend_right=True)


# ============ report mode: PDFs, faceted by sequence length ================


def tidy(df: pd.DataFrame) -> pd.DataFrame:
    """Add the display columns shared by every figure.

    ``model`` carries the typeset model name used on the x axis and ``length``
    the panel label; the raw columns are left alone so callers can still group
    on them. Both are plain strings and every plotting call passes an explicit
    order: seaborn 0.13 assigns facet data by order of first appearance rather
    than by category order, so an ordered Categorical here silently puts the
    right panel titles on the wrong panels.
    """
    df = df.copy()
    df["model"] = df["analysis_model"].map(MODEL_LABELS)
    if "sequence_length" in df.columns:
        df["length"] = df["sequence_length"].map(lambda v: f"{v} bp")
    return df


def length_order(df: pd.DataFrame) -> list:
    """Panel labels for sequence length, shortest first."""
    return [f"{v} bp" for v in sorted(df["sequence_length"].unique())]


def finalize(g: sns.FacetGrid, out_dir: Path, name: str,
             ylabel: str, xlabel: str, legend_right: bool = False,
             legend_below: int = 0, rotate_x: int = 30) -> Path:
    """Apply consistent, collision-free layout to a facet grid and save.

    Per-axes y labels are replaced by one figure-level label, x tick labels are
    rotated so model names do not run together, and an "outside" legend is
    re-anchored in figure-fraction coordinates after measuring its rendered
    width (seaborn's legend_out places it over an interior panel under
    row+col faceting).
    """
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelrotation=rotate_x)
        ax.title.set_size(12)
        ax.set_ylabel("")
        ax.set_xlabel("")
    # Margins are reserved in inches and only then converted to the figure
    # fractions subplots_adjust wants: the grids differ in size, so a fixed
    # fraction that clears the rotated tick labels on one figure crushes or
    # strands the shared labels on another.
    tick_band_in = 0.80 if rotate_x else 0.45
    bottom_in = tick_band_in + 0.35
    left_in, top_in = 1.05, 0.45
    legend_w_in = 0.15

    fig_w, fig_h = g.figure.get_size_inches()
    if legend_below and g._legend is not None:
        # A legend beside the grid has to be paid for in width, and long
        # method names make the figure wide enough that \textwidth scales it
        # down until nothing is readable. Laying it out underneath costs
        # height instead, which the page has to spare.
        # set_ncols() alone does nothing: matplotlib builds the legend's box
        # once in __init__, so the legend has to be rebuilt to re-flow it.
        sns.move_legend(g, "lower center", bbox_to_anchor=(0.5, 0.0),
                        ncol=legend_below, title=g._legend.get_title().get_text())
        g.figure.canvas.draw()
        bottom_in += g._legend.get_window_extent().height / g.figure.dpi + 0.15
    elif legend_right and g._legend is not None:
        g.figure.canvas.draw()
        legend_w_in = g._legend.get_window_extent().width / g.figure.dpi + 0.35

    fig_h += bottom_in
    g.figure.set_size_inches(fig_w + legend_w_in, fig_h)
    right_margin = fig_w / (fig_w + legend_w_in)
    if legend_below and g._legend is not None:
        g._legend.set_bbox_to_anchor((0.5, 0.0), transform=g.figure.transFigure)
        g._legend.set_loc("lower center")
    elif legend_right and g._legend is not None:
        g._legend.set_bbox_to_anchor((right_margin, 0.55), transform=g.figure.transFigure)
        g._legend.set_loc("center left")

    g.figure.supylabel(ylabel, x=0.012)
    g.figure.supxlabel(xlabel, y=(bottom_in - tick_band_in - 0.3) / fig_h
                       if legend_below else 0.012)
    g.figure.subplots_adjust(
        top=1 - top_in / fig_h, bottom=bottom_in / fig_h,
        left=left_in / (fig_w + legend_w_in), right=right_margin,
        hspace=0.45 if rotate_x else 0.3, wspace=0.22)

    # Vector PDF rather than a raster: pdfTeX copies a PDF straight through,
    # whereas an RGBA PNG (which is all matplotlib's PNG writer emits) has to
    # be decoded and split into an image plus a soft mask on every pass, which
    # alone was most of the report's compile time. CreationDate is dropped so
    # the output carries no wall-clock timestamp.
    path = out_dir / f"{name}.pdf"
    g.figure.savefig(path, metadata={"CreationDate": None})
    plt.close(g.figure)
    print(f"wrote {path}")
    return path


def report_branch_error(rep: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 1: gene-tree branch error rate (H1)."""
    df = tidy(rep)
    g = sns.catplot(
        data=df, x="model", y="branch_error_rate", hue="ils", col="length",
        kind="bar", order=model_display, hue_order=ils_order,
        col_order=length_order(df),
        height=4.0, aspect=0.95, errorbar=("ci", 95), capsize=0.15,
    )
    g.set_titles("{col_name}", size=12)
    g.legend.set_title("ILS")
    finalize(g, out_dir, "stage1_branch_error",
             ylabel="Gene-tree branch error rate", xlabel="Analysis model",
             legend_right=True)


def report_calibration_scores(rep: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 1: Brier score and expected calibration error (H2)."""
    df = tidy(rep)
    melted = df.melt(
        id_vars=["length", "ils", "model", "replicate"],
        value_vars=["brier_score", "expected_calibration_error"],
        var_name="metric", value_name="score",
    )
    metric_labels = {
        "brier_score": "Brier score",
        "expected_calibration_error": "Expected calibration error",
    }
    melted["metric"] = melted["metric"].map(metric_labels)
    g = sns.catplot(
        data=melted, x="model", y="score", hue="ils", col="length", row="metric",
        kind="bar", order=model_display, hue_order=ils_order,
        col_order=length_order(df), row_order=list(metric_labels.values()),
        height=3.6, aspect=1.1, errorbar=("ci", 95), capsize=0.15,
        sharey=False,
    )
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("ILS")
    finalize(g, out_dir, "stage1_calibration_scores",
             ylabel="Score (lower is better)", xlabel="Analysis model",
             legend_right=True)


def report_high_support_wrong(rep: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 1: how often confidently supported branches are wrong (H2)."""
    df = tidy(rep)
    melted = df.melt(
        id_vars=["length", "ils", "model", "replicate"],
        value_vars=["p_wrong_given_support_ge_0.95", "p_wrong_given_support_ge_0.99"],
        var_name="threshold", value_name="p_wrong",
    )
    threshold_labels = {
        "p_wrong_given_support_ge_0.95": "support ≥ 0.95",
        "p_wrong_given_support_ge_0.99": "support ≥ 0.99",
    }
    melted["threshold"] = melted["threshold"].map(threshold_labels)
    g = sns.catplot(
        data=melted, x="model", y="p_wrong", hue="threshold", col="length", row="ils",
        kind="bar", order=model_display, hue_order=list(threshold_labels.values()),
        col_order=length_order(df), row_order=ils_order, height=3.6, aspect=1.1,
        errorbar=("ci", 95), capsize=0.15,
    )
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("")
    finalize(g, out_dir, "stage1_high_support_error",
             ylabel="P(wrong | support ≥ t)", xlabel="Analysis model",
             legend_right=True)


def report_calibration_curve(cal: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 1: reliability curve, mean support against observed correctness."""
    df = tidy(cal)
    g = sns.relplot(
        data=df.sort_values("mean_support"), x="mean_support", y="p_correct",
        hue="model", hue_order=model_display, col="length", row="ils",
        col_order=length_order(df), row_order=ils_order,
        kind="line", marker="o", height=3.6, aspect=1.1,
        facet_kws={"sharex": True, "sharey": True},
    )
    for ax in g.axes.flat:
        ax.plot([0, 1], [0, 1], ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("Analysis model")
    finalize(g, out_dir, "stage1_calibration_curve",
             ylabel="Observed P(correct)", xlabel="Mean aBayes support in bin",
             legend_right=True, rotate_x=0)


def report_stage2_nrf(stage2: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 2: species-tree error of all four summary methods."""
    df = tidy(stage2)
    method_labels = {
        "unweighted_nrf": "ASTRAL (resolved)",
        "contracted_unweighted_nrf": "ASTRAL (0.90-contracted)",
        "weighted_nrf": "wASTRAL (support-only)",
        "hybrid_weighted_nrf": "wASTRAL (hybrid)",
    }
    melted = df.melt(
        id_vars=["length", "ils", "replicate", "model"],
        value_vars=list(method_labels), var_name="method", value_name="nrf",
    )
    melted["method"] = melted["method"].map(method_labels)
    g = sns.catplot(
        data=melted, x="model", y="nrf", hue="method", col="length", row="ils",
        kind="box", order=model_display, hue_order=list(method_labels.values()),
        col_order=length_order(df), row_order=ils_order,
        height=3.6, aspect=1.25, fliersize=2,
    )
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("Method")
    finalize(g, out_dir, "stage2_species_tree_error",
             ylabel="nRF to true species tree", xlabel="Analysis model",
             legend_below=4)


def report_stage2_deltas(stage2: pd.DataFrame, out_dir: Path, model_display: list, ils_order: list) -> None:
    """Stage 2: primary paired contrasts against the contracted baseline (H3)."""
    df = tidy(stage2)
    contrast_labels = {
        "delta_weighted_minus_contracted": "support-only wASTRAL − contracted ASTRAL",
        "delta_hybrid_minus_contracted": "hybrid wASTRAL − contracted ASTRAL",
    }
    melted = df.melt(
        id_vars=["length", "ils", "replicate", "model"],
        value_vars=list(contrast_labels), var_name="contrast", value_name="delta_nrf",
    )
    melted["contrast"] = melted["contrast"].map(contrast_labels)
    g = sns.catplot(
        data=melted, x="model", y="delta_nrf", hue="contrast", col="length", row="ils",
        kind="box", order=model_display, hue_order=list(contrast_labels.values()),
        col_order=length_order(df), row_order=ils_order,
        height=3.6, aspect=1.25, fliersize=2,
    )
    for ax in g.axes.flat:
        ax.axhline(0, ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("Contrast")
    finalize(g, out_dir, "stage2_primary_contrasts",
             ylabel="Δ nRF (negative favours wASTRAL)", xlabel="Analysis model",
             legend_below=2)


# ================================== main ===================================


def run_explore(args, models, ils, seq_lengths, designs, exclude_designs) -> None:
    """Every design on disk, faceted by design, written as PNGs."""
    stage1_summary = load_all(args.work_dir, "stage1_summary.csv", designs, exclude_designs)
    stage1_calibration = load_all(args.work_dir, "stage1_calibration.csv", designs, exclude_designs)
    stage2 = load_all(args.work_dir, "stage2_species_tree_error.csv", designs, exclude_designs)
    ils_summary = load_all(args.work_dir, "ils_summary.csv", designs, exclude_designs)

    stage1_summary = apply_filters(stage1_summary, models, ils, seq_lengths)
    stage1_calibration = apply_filters(stage1_calibration, None, ils, None)
    stage2 = apply_filters(stage2, models, ils, seq_lengths)
    ils_summary = apply_filters(ils_summary, None, ils, None)

    # These plots facet across designs; split_frames() optionally breaks each
    # into one call per design so every design gets its own PNG file.
    for suffix, (df,) in split_frames([stage1_summary], args.split_by_design):
        plot_branch_error(df, args.out, models, suffix=suffix)
        plot_calibration_metrics(df, args.out, models, suffix=suffix)
        plot_high_support_wrong(df, args.out, models, suffix=suffix)
        plot_replicate_variability(df, args.out, models, suffix=suffix)
    for suffix, (df,) in split_frames([stage1_calibration], args.split_by_design):
        plot_calibration_curve(df, args.out, models, suffix=suffix)
    for suffix, (df,) in split_frames([stage2], args.split_by_design):
        plot_stage2_nrf(df, args.out, models, suffix=suffix)
        plot_stage2_deltas(df, args.out, models, suffix=suffix)

    # These compare across designs by construction, so they are never split.
    plot_ils_separation(ils_summary, args.out)
    plot_seq_length_effect(stage1_summary, stage2, args.out, models)
    plot_model_ranking(stage1_summary, stage2, args.out, models)
    plot_pairwise_model_heatmap(stage1_summary, args.out, models)


def run_report(args, models, ils, seq_lengths, designs, exclude_designs) -> None:
    """The 50x10 designs only, faceted by sequence length, written as PDFs."""
    rep = load_all(args.work_dir, "stage1_replicate_summary.csv", designs, exclude_designs)
    cal = load_all(args.work_dir, "stage1_calibration.csv", designs, exclude_designs)
    stage2 = load_all(args.work_dir, "stage2_species_tree_error.csv", designs, exclude_designs)
    print(f"designs: {sorted(rep['design'].unique())}")

    cal = parse_calibration_condition(cal)
    rep = apply_filters(rep, models, ils, seq_lengths)
    cal = apply_filters(cal, models, ils, seq_lengths)
    stage2 = apply_filters(stage2, models, ils, seq_lengths)

    # The display orders have to match what survived filtering, otherwise
    # seaborn reserves a slot for a model or ILS level with no rows left.
    display = [MODEL_LABELS[m] for m in models if m in MODEL_LABELS]
    ils_levels = [lvl for lvl in ILS_ORDER if lvl in set(rep["ils"])]

    report_branch_error(rep, args.out, display, ils_levels)
    report_calibration_scores(rep, args.out, display, ils_levels)
    report_high_support_wrong(rep, args.out, display, ils_levels)
    report_calibration_curve(cal, args.out, display, ils_levels)
    report_stage2_nrf(stage2, args.out, display, ils_levels)
    report_stage2_deltas(stage2, args.out, display, ils_levels)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=("report", "explore"), default="report",
                        help="report: PDFs for the write-up, 50x10 designs faceted by "
                             "sequence length. explore: PNGs for every design, faceted "
                             "by design, plus the cross-design summaries. "
                             "Default: %(default)s")
    parser.add_argument("--work-dir", type=Path, default=Path("work"))
    parser.add_argument("--out", type=Path, default=None,
                        help="Output directory. Default: report/figures in report mode, "
                             "work/plots in explore mode")
    parser.add_argument("--models", type=str, default=None,
                        help="Comma-separated subset/order of models, e.g. 'gtr_g4,gtr'. "
                             f"Default: {','.join(INCLUDE_MODELS)}")
    parser.add_argument("--designs", type=str, default=None,
                        help="Comma-separated substrings; only design dirs matching one of "
                             "these are included. Default: all in explore mode, "
                             f"{','.join(REPORT_DESIGNS)} in report mode")
    parser.add_argument("--exclude-designs", type=str, default=None,
                        help="Comma-separated substrings; design dirs matching one of these "
                             "are skipped")
    parser.add_argument("--ils", type=str, default=None,
                        help="Comma-separated ILS levels to keep, e.g. 'low,high'")
    parser.add_argument("--seq-lengths", type=str, default=None,
                        help="Comma-separated sequence lengths to keep, e.g. '200,800,1600'")
    parser.add_argument("--split-by-design", action="store_true", default=SPLIT_BY_DESIGN,
                        help="Explore mode only: write one PNG per design (suffixed "
                             "__<design>) instead of faceting every design into a single "
                             "combined PNG")
    args = parser.parse_args()

    models = args.models.split(",") if args.models else INCLUDE_MODELS
    designs = args.designs.split(",") if args.designs else INCLUDE_DESIGNS
    exclude_designs = args.exclude_designs.split(",") if args.exclude_designs else EXCLUDE_DESIGNS
    ils = args.ils.split(",") if args.ils else INCLUDE_ILS
    seq_lengths = ([int(s) for s in args.seq_lengths.split(",")]
                   if args.seq_lengths else INCLUDE_SEQ_LENGTHS)
    if args.mode == "report" and designs is None:
        designs = REPORT_DESIGNS
    if args.out is None:
        args.out = DEFAULT_OUT[args.mode]

    args.out.mkdir(parents=True, exist_ok=True)

    runner = run_report if args.mode == "report" else run_explore
    runner(args, models, ils, seq_lengths, designs, exclude_designs)


if __name__ == "__main__":
    main()
