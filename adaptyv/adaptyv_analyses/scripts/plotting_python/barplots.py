#!/usr/bin/env python3
"""
Standalone Matplotlib barplots script matching R barplots.R behavior.

Usage:
  python scripts/plotting_python/barplots.py \
    --input ./data/processed/all_submissions_new.csv \
    --x_column design_category \
    --color_column selected \
    --output ./plots/barplots \
    --format svg --width 2600 --height 2200 --dpi 300 \
    --round both --title "Design category" \
    --top_n 10 --subtitle "Number of designs submitted vs. selected for validation" \
    --sort size
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # type: ignore


def _ensure_roboto_installed() -> str:
  try:
    from matplotlib import font_manager
    available = {f.name for f in font_manager.fontManager.ttflist}
    if "Roboto" in available:
      return "Roboto"
  except Exception:
    pass
  return "DejaVu Sans"


def set_theme() -> None:
  roboto = _ensure_roboto_installed()
  mpl.rcParams.update({
    "font.family": roboto,
    "text.color": "#333333",
    "axes.labelcolor": "#333333",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": False,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.5,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.framealpha": 1.0,
    "legend.facecolor": "white",
    "axes.titlesize": 20,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
  })


def get_palettes() -> Dict[str, Dict[str, str]]:
  binding_colors = {
    "None": "#E5E7EB",
    "Not expressed": "#9EA2AF",
    "Weak": "#C4B2FB",
    "Medium": "#A7C1FB",
    "Strong": "#2A4DD0",
    "Missing binding data": "#3E6175",
    "% binders": "#DC7A73",
  }
  expression_colors = {"None": "#E5E7EB", "Low": "#9EA2AF", "Medium": "#E9C435", "High": "#69B7A7"}
  selection_status_colors = {"Top 100": "#8B90DD", "Adaptyv selection": "#8CD2F4", "Not selected": "#E5E7EB"}
  design_category_colors = {"De novo": "#56A6D4", "Optimized binder": "#8B90DD", "Diversified binder": "#8CD2F4", "Hallucination": "#3E6175"}
  bindcraft_colors = {"BindCraft": "#56A6D4", "Other": "#E5E7EB"}
  rfdiffusion_colors = {"RFdiffusion": "#8B90DD", "Other": "#E5E7EB"}
  esm_colors = {"ESM": "#8CD2F4", "Other": "#E5E7EB"}
  timed_colors = {"TIMED": "#56A6D4", "ProteinMPNN": "#8B90DD", "ESM-IF": "#8CD2F4", "Other": "#E5E7EB"}
  af2_backprop_colors = {"AF2 backprop": "#8CD2F4", "Other": "#E5E7EB"}
  other_hallucination_colors = {"Other hallucination": "#8CD2F4", "Other": "#E5E7EB"}
  metric_colors = {"ESM2 PLL": "#8CD2F4", "ipTM": "#8B90DD", "iPAE": "#3E6175"}

  all_colors_list: List[str] = []
  for pal in [binding_colors, expression_colors, selection_status_colors, design_category_colors,
              bindcraft_colors, rfdiffusion_colors, esm_colors, timed_colors,
              af2_backprop_colors, other_hallucination_colors, metric_colors]:
    all_colors_list.extend(list(pal.values()))
  all_unique = list(dict.fromkeys(all_colors_list))

  return {
    "binding_strength": binding_colors,
    "expression": expression_colors,
    "selected": selection_status_colors,
    "design_category": design_category_colors,
    "BindCraft": bindcraft_colors,
    "RFdiffusion": rfdiffusion_colors,
    "ESM": esm_colors,
    "TIMED": timed_colors,
    "AF2_backprop": af2_backprop_colors,
    "Other_hallucination": other_hallucination_colors,
    "metric": metric_colors,
    "_all": {str(i): c for i, c in enumerate(all_unique)},
  }


def apply_blog_post_theme(fig: plt.Figure, ax: plt.Axes, *, title: str, subtitle: str, x_label: str, y_label: str, legend_title: str) -> None:
  set_theme()
  ax.set_xlabel(x_label)
  ax.set_ylabel(y_label)
  if title:
    fig.suptitle(title, fontsize=20, fontweight="bold", y=0.98)
  if subtitle:
    ax.set_title(subtitle, color="#3D7A9A", fontsize=16, pad=14)
  leg = ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
  leg.get_frame().set_edgecolor("black")
  leg.get_frame().set_facecolor("white")
  plt.subplots_adjust(right=0.78, left=0.12, top=0.90, bottom=0.22)


def plot_bar(
  df: pd.DataFrame,
  *,
  x_column: str,
  color_column: str,
  output_dir: Optional[Path],
  output_format: str,
  width: int,
  height: int,
  dpi: int,
  remove_not_mentioned: bool,
  remove_missing_data: bool,
  round_filter: str,
  title: str,
  top_n: int,
  subtitle: str,
  sort: str,
) -> Tuple[plt.Figure, plt.Axes]:
  set_theme()
  data = df.copy()

  if round_filter != "both" and "round" in data.columns:
    try:
      rn = int(round_filter)
      data = data[data["round"] == rn]
    except Exception:
      pass

  if color_column in data.columns:
    data[color_column] = data[color_column].where(~data[color_column].isna(), "missing data")

  if top_n and top_n > 0 and x_column in data.columns:
    counts_by_x = data[x_column].value_counts(dropna=False)
    top_x_vals = list(counts_by_x.index[: min(top_n, counts_by_x.size)])
    data = data[data[x_column].isin(top_x_vals)]

  if remove_not_mentioned:
    if color_column in data.columns:
      data = data[data[color_column] != "Not mentioned"]
    if x_column in data.columns:
      data = data[data[x_column] != "Not mentioned"]

  if remove_missing_data:
    if color_column in data.columns:
      mask = data[color_column].notna() & (data[color_column].astype(str) != "") & (data[color_column] != "Missing data")
      data = data[mask]
    if x_column in data.columns:
      mask_x = data[x_column].notna() & (data[x_column].astype(str) != "") & (data[x_column] != "Missing data")
      data = data[mask_x]

  if x_column in data.columns and x_column == "selected":
    data = data[data[x_column] != "No"]
  if x_column in data.columns and x_column == "binding":
    data = data[data[x_column] != "Unknown"]

  palettes = get_palettes()
  if color_column == "binding_strength" and color_column in data.columns:
    data.loc[data[color_column] == "Unknown", color_column] = "Not expressed"
    color_levels = ["None", "Not expressed", "Weak", "Medium", "Strong"]
  elif color_column == "expression" and color_column in data.columns:
    color_levels = ["None", "Low", "Medium", "High"]
  elif color_column == "selected" and color_column in data.columns:
    data.loc[data[color_column] == "No", color_column] = "Not selected"
    color_levels = ["Top 100", "Adaptyv selection", "Not selected"]
  elif color_column == "BindCraft" and color_column in data.columns:
    data[color_column] = data[color_column].map({"Yes": "BindCraft", "No": "Other"}).fillna(data[color_column])
    color_levels = ["BindCraft", "Other"]
  elif color_column == "RFdiffusion" and color_column in data.columns:
    data[color_column] = data[color_column].map({"Yes": "RFdiffusion", "No": "Other"}).fillna(data[color_column])
    color_levels = ["RFdiffusion", "Other"]
  elif color_column == "ESM" and color_column in data.columns:
    data[color_column] = data[color_column].map({"Yes": "ESM", "No": "Other"}).fillna(data[color_column])
    color_levels = ["ESM", "Other"]
  elif color_column == "TIMED" and color_column in data.columns:
    valid_categories = ["TIMED", "ProteinMPNN", "ESM-IF", "Other"]
    data.loc[~data[color_column].isin(valid_categories), color_column] = "Other"
    color_levels = valid_categories
  elif color_column == "AF2_backprop" and color_column in data.columns:
    data[color_column] = data[color_column].map({"Yes": "AF2 backprop", "No": "Other"}).fillna(data[color_column])
    color_levels = ["AF2 backprop", "Other"]
  elif color_column == "Other_hallucination" and color_column in data.columns:
    data[color_column] = data[color_column].map({"Yes": "Other hallucination", "No": "Other"}).fillna(data[color_column])
    color_levels = ["Other hallucination", "Other"]
  else:
    color_levels = list(pd.unique(data[color_column])) if color_column in data.columns else []

  if sort == "size" and x_column in data.columns:
    x_order = data.groupby(x_column).size().sort_values(ascending=False).index.tolist()
  elif x_column == "binding":
    x_order = ["Yes", "No"]
  else:
    x_order = sorted(list(pd.unique(data[x_column]))) if x_column in data.columns else []

  grouped = data.groupby([x_column, color_column]).size().reset_index(name="count")
  if x_order:
    grouped[x_column] = pd.Categorical(grouped[x_column], categories=x_order, ordered=True)
  if color_levels:
    grouped[color_column] = pd.Categorical(grouped[color_column], categories=color_levels, ordered=True)
  grouped = grouped.sort_values([x_column, color_column])

  totals = grouped.groupby(x_column, observed=False)["count"].sum().reset_index(name="total")

  def _percent_for_row(x_value: str) -> int:
    sub = grouped[grouped[x_column] == x_value]
    total_n = int(sub["count"].sum())
    if total_n == 0:
      return 0
    if color_column == "selected":
      sel = sub[sub[color_column].isin(["Top 100", "Adaptyv selection"])]
      return int(round(100 * sel["count"].sum() / total_n))
    if color_column == "binding_strength":
      sel = sub[sub[color_column].isin(["Weak", "Medium", "Strong"])]
      return int(round(100 * sel["count"].sum() / total_n))
    if color_column == "expression":
      sel = sub[sub[color_column].isin(["Low", "Medium", "High"])]
      return int(round(100 * sel["count"].sum() / total_n))
    if color_column == "TIMED":
      sel = sub[sub[color_column] != "Other"]
      return int(round(100 * sel["count"].sum() / total_n))
    if color_column in {"BindCraft", "RFdiffusion", "ESM", "AF2_backprop", "Other_hallucination"}:
      sel = sub[sub[color_column].isin(["BindCraft", "RFdiffusion", "ESM", "AF2 backprop", "Other hallucination"])]
      non_other = sel[~sel[color_column].isin(["Other"])]
      return int(round(100 * non_other["count"].sum() / total_n))
    return 0

  totals["selected"] = [_percent_for_row(val) for val in totals[x_column].astype(str).tolist()]

  palette_all = get_palettes()
  palette = palette_all.get(color_column, palette_all["_all"])
  categories_present = [c for c in color_levels if c in grouped[color_column].unique().tolist()]
  color_map: Dict[str, str] = {}
  fallback_colors = list(palette_all["_all"].values())
  idx = 0
  for cat in categories_present:
    color_map[cat] = palette.get(cat, fallback_colors[idx % len(fallback_colors)])
    idx += 1

  fig_w_inches = max(1.0, width / float(dpi))
  fig_h_inches = max(1.0, height / float(dpi))
  fig, ax = plt.subplots(figsize=(fig_w_inches, fig_h_inches), dpi=dpi)

  bar_x = np.arange(len(x_order)) if x_order else np.arange(len(totals))
  bottoms = np.zeros(len(bar_x), dtype=float)

  for level in color_levels:
    if level not in grouped[color_column].cat.categories:
      continue
    sub = grouped[grouped[color_column] == level]
    counts: List[int] = []
    for xv in x_order:
      row = sub[sub[x_column] == xv]
      counts.append(int(row["count"].iloc[0]) if not row.empty else 0)
    counts_arr = np.array(counts, dtype=float)
    bars = ax.bar(bar_x, counts_arr, width=0.8, bottom=bottoms, color=color_map.get(level, "#E5E7EB"), edgecolor="black", linewidth=1.0, label=str(level))
    for rect, c in zip(bars, counts_arr):
      if c > 0 and rect.get_height() >= 1:
        ax.text(rect.get_x() + rect.get_width()/2.0, rect.get_y() + rect.get_height()/2.0, f"{int(c)}", ha="center", va="center", color="white", fontsize=10, clip_on=True)
    bottoms += counts_arr

  totals = totals.set_index(x_column).loc[x_order].reset_index()
  percent_suffix = (
    "binders" if color_column == "binding_strength" else
    "expressed" if color_column == "expression" else
    "selected" if color_column == "selected" else
    "design methods" if color_column == "TIMED" else
    color_column if color_column in {"BindCraft", "RFdiffusion", "ESM", "AF2_backprop", "Other_hallucination"} else ""
  )
  for idx2, (_, total_val, sel_pct) in enumerate(totals[[x_column, "total", "selected"]].itertuples(index=False, name=None)):
    label_txt = f"{sel_pct}% {percent_suffix}" if percent_suffix else f"{sel_pct}%"
    ax.text(bar_x[idx2], float(total_val) + max(0.02 * float(total_val), 0.1), label_txt, ha="center", va="bottom", color="#3D7A9A", fontsize=12, fontweight="bold")

  ax.set_xticks(bar_x)
  ax.set_xticklabels([str(x) for x in x_order], rotation=45, ha="right")

  legend_title = ("Design model" if color_column in {"BindCraft", "RFdiffusion", "ESM", "TIMED", "AF2_backprop", "Other_hallucination"}
                  else x_column.replace("_", " ").title())
  apply_blog_post_theme(fig, ax, title=title, subtitle=subtitle, x_label=x_column.replace("_", " ").title(), y_label="Number of designs", legend_title=legend_title)

  ymin, ymax = ax.get_ylim()
  ax.set_ylim(ymin, ymax * 1.1 if ymax > 0 else 1)

  if output_dir is not None:
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = output_dir / f"barplot_{x_column}_by_{color_column}.{output_format}"
    if output_format.lower() == "png":
      fig.savefig(filename, dpi=dpi, facecolor=fig.get_facecolor(), edgecolor="none")
    elif output_format.lower() == "svg":
      fig.savefig(filename, format="svg", facecolor="none")
    else:
      fig.savefig(filename)

  return fig, ax


def main() -> None:
  p = argparse.ArgumentParser()
  p.add_argument("--input", type=Path, default=Path("./data/processed/all_submissions_new.csv"))
  p.add_argument("--x_column", type=str, default="design_category")
  p.add_argument("--color_column", type=str, default="selected")
  p.add_argument("--output", type=Path, default=Path("./plots/barplots"))
  p.add_argument("--format", type=str, default="svg")
  p.add_argument("--width", type=int, default=2600)
  p.add_argument("--height", type=int, default=2200)
  p.add_argument("--dpi", type=int, default=300)
  p.add_argument("--remove_not_mentioned", action="store_true")
  p.add_argument("--remove_missing_data", action="store_true")
  p.add_argument("--round", dest="round_filter", type=str, default="both")
  p.add_argument("--title", type=str, default="Design category")
  p.add_argument("--top_n", type=int, default=10)
  p.add_argument("--subtitle", type=str, default="Number of designs submitted vs. selected for validation")
  p.add_argument("--sort", type=str, default="size", choices=["size", "alpha"])
  args = p.parse_args()

  df = pd.read_csv(args.input)
  plot_bar(
    df,
    x_column=args.x_column,
    color_column=args.color_column,
    output_dir=args.output,
    output_format=args.format,
    width=args.width,
    height=args.height,
    dpi=args.dpi,
    remove_not_mentioned=args.remove_not_mentioned,
    remove_missing_data=args.remove_missing_data,
    round_filter=args.round_filter,
    title=args.title,
    top_n=args.top_n,
    subtitle=args.subtitle,
    sort=args.sort,
  )


if __name__ == "__main__":
  main()


