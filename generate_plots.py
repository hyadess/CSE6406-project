#!/usr/bin/env python3
"""Generate report figures from the 50x10 experiments under work/*/results.

Usage:
    .venv-plots/bin/python generate_plots.py [--work-dir work] [--out report/figures]

Only the ``paper_subset_50x10_*`` design directories are read: these are the
full-size runs (50 loci, 10 replicates) and they partition cleanly by sequence
length, one design per length. Rows are therefore tagged with their sequence
length and every figure facets on that rather than on the design name.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DESIGN_GLOB = "paper_subset_50x10_*"

MODEL_ORDER = ["gtr_g4", "gtr", "hky_g4", "hky", "jc"]
MODEL_LABELS = {
    "gtr_g4": "GTR+Γ4",
    "gtr": "GTR",
    "hky_g4": "HKY+Γ4",
    "hky": "HKY",
    "jc": "JC",
}
MODEL_DISPLAY = [MODEL_LABELS[m] for m in MODEL_ORDER]
ILS_ORDER = ["low", "high"]

sns.set_theme(style="whitegrid", context="talk")
matplotlib.rcParams["pdf.fonttype"] = 42


def load_all(work_dir: Path, csv_name: str, design_glob: str) -> pd.DataFrame:
    """Concatenate one CSV across every matching design directory."""
    frames = []
    for design_dir in sorted(work_dir.glob(design_glob)):
        csv_path = design_dir / "results" / csv_name
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        df.insert(0, "design", design_dir.name)
        frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"No {csv_name} found under {work_dir}/{design_glob}/results")
    return pd.concat(frames, ignore_index=True)


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


def plot_branch_error(rep: pd.DataFrame, out_dir: Path) -> None:
    """Stage 1: gene-tree branch error rate (H1)."""
    df = tidy(rep)
    g = sns.catplot(
        data=df, x="model", y="branch_error_rate", hue="ils", col="length",
        kind="bar", order=MODEL_DISPLAY, hue_order=ILS_ORDER,
        col_order=length_order(df),
        height=4.0, aspect=0.95, errorbar=("ci", 95), capsize=0.15,
    )
    g.set_titles("{col_name}", size=12)
    g.legend.set_title("ILS")
    finalize(g, out_dir, "stage1_branch_error",
             ylabel="Gene-tree branch error rate", xlabel="Analysis model",
             legend_right=True)


def plot_calibration_scores(rep: pd.DataFrame, out_dir: Path) -> None:
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
        kind="bar", order=MODEL_DISPLAY, hue_order=ILS_ORDER,
        col_order=length_order(df), row_order=list(metric_labels.values()),
        height=3.6, aspect=1.1, errorbar=("ci", 95), capsize=0.15,
        sharey=False,
    )
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("ILS")
    finalize(g, out_dir, "stage1_calibration_scores",
             ylabel="Score (lower is better)", xlabel="Analysis model",
             legend_right=True)


def plot_high_support_wrong(rep: pd.DataFrame, out_dir: Path) -> None:
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
        kind="bar", order=MODEL_DISPLAY, hue_order=list(threshold_labels.values()),
        col_order=length_order(df), row_order=ILS_ORDER, height=3.6, aspect=1.1,
        errorbar=("ci", 95), capsize=0.15,
    )
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("")
    finalize(g, out_dir, "stage1_high_support_error",
             ylabel="P(wrong | support ≥ t)", xlabel="Analysis model",
             legend_right=True)


def plot_calibration_curve(cal: pd.DataFrame, out_dir: Path) -> None:
    """Stage 1: reliability curve, mean support against observed correctness."""
    df = cal.copy()
    df["analysis_model"] = df["condition"].str.extract(r"model=([a-z_0-9]+)")
    df["ils"] = df["condition"].str.extract(r"ils=([a-z]+)")
    df["sequence_length"] = df["condition"].str.extract(r"length=(\d+)").astype(int)
    df = tidy(df)
    g = sns.relplot(
        data=df.sort_values("mean_support"), x="mean_support", y="p_correct",
        hue="model", hue_order=MODEL_DISPLAY, col="length", row="ils",
        col_order=length_order(df), row_order=ILS_ORDER,
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


def plot_stage2_nrf(stage2: pd.DataFrame, out_dir: Path) -> None:
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
        kind="box", order=MODEL_DISPLAY, hue_order=list(method_labels.values()),
        col_order=length_order(df), row_order=ILS_ORDER,
        height=3.6, aspect=1.25, fliersize=2,
    )
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("Method")
    finalize(g, out_dir, "stage2_species_tree_error",
             ylabel="nRF to true species tree", xlabel="Analysis model",
             legend_below=4)


def plot_stage2_deltas(stage2: pd.DataFrame, out_dir: Path) -> None:
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
        kind="box", order=MODEL_DISPLAY, hue_order=list(contrast_labels.values()),
        col_order=length_order(df), row_order=ILS_ORDER,
        height=3.6, aspect=1.25, fliersize=2,
    )
    for ax in g.axes.flat:
        ax.axhline(0, ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="ILS = {row_name}", col_template="{col_name}", size=12)
    g.legend.set_title("Contrast")
    finalize(g, out_dir, "stage2_primary_contrasts",
             ylabel="Δ nRF (negative favours wASTRAL)", xlabel="Analysis model",
             legend_below=2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path("work"))
    parser.add_argument("--out", type=Path, default=Path("report/figures"))
    parser.add_argument("--design-glob", default=DESIGN_GLOB,
                        help="which design directories to include (default: %(default)s)")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    rep = load_all(args.work_dir, "stage1_replicate_summary.csv", args.design_glob)
    cal = load_all(args.work_dir, "stage1_calibration.csv", args.design_glob)
    stage2 = load_all(args.work_dir, "stage2_species_tree_error.csv", args.design_glob)
    print(f"designs: {sorted(rep['design'].unique())}")

    plot_branch_error(rep, args.out)
    plot_calibration_scores(rep, args.out)
    plot_high_support_wrong(rep, args.out)
    plot_calibration_curve(cal, args.out)
    plot_stage2_nrf(stage2, args.out)
    plot_stage2_deltas(stage2, args.out)


if __name__ == "__main__":
    main()
