# pyright: reportMissingTypeStubs=false
"""
Matplotlib theme and style utilities matching the R ggplot theme.

Includes:
- set_adaptyv_matplotlib_theme(): global rcParams to match R theme
- get_adaptyv_palettes(): color palettes used across plots
- apply_adaptyv_blog_post_theme(): apply titles/legend/layout to an Axes
"""

from __future__ import annotations

from typing import Dict, List, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


# -----------------------------
# Theme configuration
# -----------------------------

def _ensure_roboto_installed() -> str:
  """Return the name of a usable Roboto font family if available, else fallback.

  Attempts to use "Roboto" if present in the Matplotlib font manager. If not
  available, returns a generic sans-serif family. This keeps plotting robust
  without requiring a runtime font download.
  """
  try:
    from matplotlib import font_manager
    available = {f.name for f in font_manager.fontManager.ttflist}
    if "Roboto" in available:
      return "Roboto"
  except Exception:
    pass
  return "DejaVu Sans"


def set_adaptyv_matplotlib_theme() -> None:
  """Apply Matplotlib rcParams to match the R `adaptyv_theme` as closely as possible.

  Mirrors scripts/plotting_r/barplots.R theme:
  - Title size 20 bold, subtitle size 16 colored #3D7A9A
  - Axis titles size 12 bold; tick labels size 10; x tick labels typically 14 in plots
  - White backgrounds; no panel grid; black axis spines 0.5 linewidth
  - Legend on right with white background and black border
  """
  roboto = _ensure_roboto_installed()

  mpl.rcParams.update({
    # Base font
    "font.family": roboto,
    "text.color": "#333333",
    "axes.labelcolor": "#333333",
    "xtick.color": "#333333",
    "ytick.color": "#333333",

    # Figure and axes backgrounds
    "figure.facecolor": "white",
    "axes.facecolor": "white",

    # Grid (off, to match panel.grid.blank in R barplots.R)
    "axes.grid": False,

    # Spines (axis lines)
    "axes.edgecolor": "black",
    "axes.linewidth": 0.5,

    # Legend styling
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.framealpha": 1.0,
    "legend.facecolor": "white",

    # Title/labels sizing defaults (exact values set in plotting function)
    "axes.titlesize": 20,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
  })


# -----------------------------
# Palettes (mirroring barplots.R)
# -----------------------------

def get_adaptyv_palettes() -> Dict[str, Dict[str, str]]:
  """Return all palettes used across plots, matching R definitions exactly."""
  binding_colors = {
    "None": "#E5E7EB",
    "Not expressed": "#9EA2AF",
    "Weak": "#C4B2FB",
    "Medium": "#A7C1FB",
    "Strong": "#2A4DD0",
    "Missing binding data": "#3E6175",
    "% binders": "#DC7A73",
  }

  expression_colors = {
    "None": "#E5E7EB",
    "Low": "#9EA2AF",
    "Medium": "#E9C435",
    "High": "#69B7A7",
  }

  selection_status_colors = {
    "Top 100": "#8B90DD",
    "Adaptyv selection": "#8CD2F4",
    "Not selected": "#E5E7EB",
  }

  design_category_colors = {
    "De novo": "#56A6D4",
    "Optimized binder": "#8B90DD",
    "Diversified binder": "#8CD2F4",
    "Hallucination": "#3E6175",
  }

  # Binary model-type columns map Yes/No to distinctive vs. gray
  model_type_colors = {"Yes": "#56A6D4", "No": "#E5E7EB"}

  bindcraft_colors = {"BindCraft": "#56A6D4", "Other": "#E5E7EB"}
  rfdiffusion_colors = {"RFdiffusion": "#8B90DD", "Other": "#E5E7EB"}
  esm_colors = {"ESM": "#8CD2F4", "Other": "#E5E7EB"}

  timed_colors = {
    "TIMED": "#56A6D4",
    "ProteinMPNN": "#8B90DD",
    "ESM-IF": "#8CD2F4",
    "Other": "#E5E7EB",
  }

  af2_backprop_colors = {"AF2 backprop": "#8CD2F4", "Other": "#E5E7EB"}
  other_hallucination_colors = {"Other hallucination": "#8CD2F4", "Other": "#E5E7EB"}

  metric_colors = {"ESM2 PLL": "#8CD2F4", "ipTM": "#8B90DD", "iPAE": "#3E6175"}

  # Aggregate unique palette colors for fallback
  all_colors_list: List[str] = []
  for pal in [
    binding_colors,
    expression_colors,
    selection_status_colors,
    design_category_colors,
    model_type_colors,
    bindcraft_colors,
    rfdiffusion_colors,
    esm_colors,
    timed_colors,
    af2_backprop_colors,
    other_hallucination_colors,
    metric_colors,
  ]:
    all_colors_list.extend(list(pal.values()))
  adaptyv_colors_unique = list(dict.fromkeys(all_colors_list))

  return {
    "binding_strength": binding_colors,
    "expression": expression_colors,
    "selected": selection_status_colors,
    "design_category": design_category_colors,
    "model_type": model_type_colors,
    "BindCraft": bindcraft_colors,
    "RFdiffusion": rfdiffusion_colors,
    "ESM": esm_colors,
    "TIMED": timed_colors,
    "AF2_backprop": af2_backprop_colors,
    "Other_hallucination": other_hallucination_colors,
    "metric": metric_colors,
    "_all": {str(i): c for i, c in enumerate(adaptyv_colors_unique)},
  }


def _palette_for_column(color_column: str, palettes: Dict[str, Dict[str, str]]) -> Dict[str, str]:
  """Return palette for a given column, or the combined palette as fallback."""
  return palettes.get(color_column, palettes["_all"])


# -----------------------------
# Blog post theme application (titles, legend, layout)
# -----------------------------

def apply_adaptyv_blog_post_theme(
  fig: Figure,
  ax: Axes,
  *,
  title: str,
  subtitle: str,
  x_label: Optional[str] = None,
  y_label: Optional[str] = "Number of designs",
  legend_title: Optional[str] = None,
) -> None:
  """Apply blog-post styling to the axes and figure.

  - Main title as figure suptitle to avoid overlap
  - Subtitle as axes title in #3D7A9A
  - Legend on the right with black border and white bg
  - Adjust layout to make room for legend and tick labels
  """
  set_adaptyv_matplotlib_theme()

  if x_label is not None:
    ax.set_xlabel(x_label)
  if y_label is not None:
    ax.set_ylabel(y_label)

  if title:
    fig.suptitle(title, fontsize=20, fontweight="bold", y=0.98)
  if subtitle:
    ax.set_title(subtitle, color="#3D7A9A", fontsize=16, pad=14)

  if legend_title is not None:
    leg = ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
  else:
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
  if leg is not None:
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_facecolor("white")

  plt.subplots_adjust(right=0.78, left=0.12, top=0.90, bottom=0.22)


# Plotting functions have been moved to scripts/plotting_python/barplots.py.


__all__ = [
  "set_adaptyv_matplotlib_theme",
  "apply_adaptyv_blog_post_theme",
  "get_adaptyv_palettes",
]

 
