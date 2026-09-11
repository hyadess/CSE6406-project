#!/usr/bin/env python3
"""Generate comparison figures from the CSVs under work/*/results.

Usage:
    .venv-plots/bin/python generate_plots.py [--work-dir work] [--out plots]

Reads every ``paper_subset_*`` (or any) design directory under --work-dir that
contains a ``results/`` folder, tags rows with the design name, concatenates
across designs, and writes a set of PNG figures comparing analysis models,
ILS levels, and Stage 2 species-tree methods.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

MODEL_ORDER = ["gtr_g4", "gtr", "hky_g4", "hky", "jc"]
sns.set_theme(style="whitegrid", context="talk")


def load_all(work_dir: Path, csv_name: str) -> pd.DataFrame:
    frames = []
    for design_dir in sorted(work_dir.iterdir()):
        csv_path = design_dir / "results" / csv_name
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        df.insert(0, "design", design_dir.name)
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No {csv_name} found under {work_dir}/*/results")
    return pd.concat(frames, ignore_index=True)


def order_models(df: pd.DataFrame, col: str = "analysis_model") -> pd.DataFrame:
    df[col] = pd.Categorical(df[col], categories=MODEL_ORDER, ordered=True)
    return df.sort_values(col)


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


def plot_branch_error(df: pd.DataFrame, out_dir: Path) -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy())
    g = sns.catplot(
        data=pooled, x="analysis_model", y="branch_error_rate", hue="ils",
        col="design", kind="bar", height=4.5, aspect=0.9, col_wrap=2,
        order=MODEL_ORDER, legend=True,
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
    g.figure.savefig(out_dir / "stage1_branch_error_by_model.png", dpi=150)
    plt.close(g.figure)
    print(f"wrote {out_dir / 'stage1_branch_error_by_model.png'}")


def plot_calibration_metrics(df: pd.DataFrame, out_dir: Path) -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy())
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sns.barplot(data=pooled, x="analysis_model", y="brier_score", hue="ils",
                order=MODEL_ORDER, ax=axes[0])
    axes[0].set_title("Brier score (lower = better)")
    axes[0].set_xlabel("Analysis model")
    sns.barplot(data=pooled, x="analysis_model", y="expected_calibration_error", hue="ils",
                order=MODEL_ORDER, ax=axes[1])
    axes[1].set_title("Expected calibration error")
    axes[1].set_xlabel("Analysis model")
    save(fig, out_dir, "stage1_calibration_scores")


def plot_high_support_wrong(df: pd.DataFrame, out_dir: Path) -> None:
    pooled = df[df["replicate"] == "pooled"] if df["replicate"].dtype == object else df
    pooled = order_models(pooled.copy())
    melted = pooled.melt(
        id_vars=["design", "ils", "analysis_model"],
        value_vars=["p_wrong_given_support_ge_0.95", "p_wrong_given_support_ge_0.99"],
        var_name="threshold", value_name="p_wrong",
    )
    g = sns.catplot(
        data=melted, x="analysis_model", y="p_wrong", hue="threshold",
        col="ils", row="design", kind="bar", order=MODEL_ORDER, height=3.5, aspect=1.2,
        legend_out=True,
    )
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, "stage1_high_support_error",
                        ylabel="P(wrong | support >= t)", xlabel="Analysis model",
                        legend_right=True)


def plot_calibration_curve(df: pd.DataFrame, out_dir: Path) -> None:
    df = df.copy()
    df["model"] = df["condition"].str.extract(r"model=([a-z_0-9]+)")
    df["ils"] = df["condition"].str.extract(r"ils=([a-z]+)")
    df["model"] = pd.Categorical(df["model"], categories=MODEL_ORDER, ordered=True)
    g = sns.relplot(
        data=df.sort_values("mean_support"), x="mean_support", y="p_correct",
        hue="model", hue_order=MODEL_ORDER, col="ils", row="design",
        kind="line", marker="o", height=3.5, aspect=1.2, facet_kws={"sharex": True, "sharey": True},
    )
    for ax in g.axes.flat:
        ax.plot([0, 1], [0, 1], ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, "stage1_calibration_curve",
                        ylabel="Observed P(correct)", xlabel="Mean support in bin",
                        legend_right=True)


def plot_stage2_nrf(df: pd.DataFrame, out_dir: Path) -> None:
    df = order_models(df.copy())
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
    g = sns.catplot(
        data=melted, x="analysis_model", y="nrf", hue="method",
        col="ils", row="design", kind="box", order=MODEL_ORDER, height=3.5, aspect=1.3,
        legend_out=True,
    )
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, "stage2_species_tree_error",
                        ylabel="Normalized RF to true species tree", xlabel="Analysis model",
                        legend_right=True)


def plot_stage2_deltas(df: pd.DataFrame, out_dir: Path) -> None:
    df = order_models(df.copy())
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
    g = sns.catplot(
        data=melted, x="analysis_model", y="delta_nrf", hue="contrast",
        col="ils", row="design", kind="box", order=MODEL_ORDER, height=3.5, aspect=1.3,
        legend_out=True,
    )
    for ax in g.axes.flat:
        ax.axhline(0, ls="--", color="grey", linewidth=1)
    g.set_titles(row_template="{row_name}", col_template="ILS = {col_name}", size=11)
    finalize_facetgrid(g, out_dir, "stage2_primary_contrasts",
                        ylabel="Delta nRF (negative favors wASTRAL)", xlabel="Analysis model",
                        legend_right=True)


def plot_ils_separation(df: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sns.barplot(data=df, x="design", y="mean_true_gene_species_nrf", hue="ils", ax=ax)
    ax.set_ylabel("Mean true gene/species nRF")
    ax.set_xlabel("Design")
    ax.tick_params(axis="x", rotation=20)
    save(fig, out_dir, "ils_separation_check")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path("work"))
    parser.add_argument("--out", type=Path, default=Path("work/plots"))
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    stage1_summary = load_all(args.work_dir, "stage1_summary.csv")
    stage1_calibration = load_all(args.work_dir, "stage1_calibration.csv")
    stage2 = load_all(args.work_dir, "stage2_species_tree_error.csv")
    ils_summary = load_all(args.work_dir, "ils_summary.csv")

    plot_branch_error(stage1_summary, args.out)
    plot_calibration_metrics(stage1_summary, args.out)
    plot_high_support_wrong(stage1_summary, args.out)
    plot_calibration_curve(stage1_calibration, args.out)
    plot_stage2_nrf(stage2, args.out)
    plot_stage2_deltas(stage2, args.out)
    plot_ils_separation(ils_summary, args.out)


if __name__ == "__main__":
    main()