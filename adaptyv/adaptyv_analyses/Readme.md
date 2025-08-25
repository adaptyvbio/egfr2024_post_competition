Adaptyv Competition – Annotation and plotting

Directory layout
- src: core code for binding property analysis and utilities
- scripts/plotting_r: R analysis and plotting scripts
- scripts/plotting_python: Python plotting utilities (Matplotlib)
  - barplots.py: standalone CLI for stacked barplots
  - blog_post_theme.py: theme and palettes (no plotting)
- data: raw, processed, and intermediate datasets
- plots: generated figures

Python setup
1) Create environment and install deps
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

2) Generate barplots (Python CLI)
```bash
python scripts/plotting_python/barplots.py \
  --input ./data/processed/all_submissions_new.csv \
  --x_column design_category \
  --color_column selected \
  --output ./plots/barplots \
  --format svg --width 2600 --height 2200 --dpi 300 \
  --round both --title "Design category" \
  --top_n 10 --subtitle "Number of designs submitted vs. selected for validation" \
  --sort size
```

Optional programmatic use of the theme/palettes
```python
from scripts.plotting_python.blog_post_theme import (
  set_adaptyv_matplotlib_theme,
  get_adaptyv_palettes,
  apply_adaptyv_blog_post_theme,
)
import matplotlib.pyplot as plt

set_adaptyv_matplotlib_theme()
palettes = get_adaptyv_palettes()
fig, ax = plt.subplots()
# ... draw your plot ...
apply_adaptyv_blog_post_theme(
  fig, ax, title="My Title", subtitle="My Subtitle",
  x_label="X", y_label="Y", legend_title="Legend"
)
fig.savefig("./plots/example.svg")
```

R setup
Install R packages used in the scripts:
```bash
Rscript ./install_packages.R
```

Run R barplots
```bash
Rscript scripts/plotting_r/barplots.R \
  --input ./data/processed/all_submissions_new.csv \
  --x_column design_category \
  --color_column selected \
  --output ./plots/barplots \
  --format svg \
  --width 2600 --height 2200 --res 300 \
  --round both --title "Design category" \
  --top_n 10 --subtitle "Number of designs submitted vs. selected for validation" \
  --sort size
```

Optional programmatic use of the theme/palettes (R)
```r
# Load theme and palettes
source("scripts/plotting_r/blog_post_theme.R")
library(ggplot2)

# Example plot using the theme and a palette
df <- data.frame(
  x = factor(c("A","B","C","D"), levels = c("A","B","C","D")),
  y = c(10, 15, 7, 12),
  group = c("De novo", "Optimized binder", "Diversified binder", "Hallucination")
)

ggplot(df, aes(x = x, y = y, fill = group)) +
  geom_col(color = "black") +
  scale_fill_manual(values = design_category_colors, name = "Design category") +
  labs(
    title = "My Title",
    subtitle = "My Subtitle",
    x = "X",
    y = "Y"
  ) +
  adaptyv_theme()
```


Amino acid composition (R)
```bash
Rscript scripts/plotting_r/aa_composition.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/aa_composition \
  --format svg --width 8400 --height 3200 --res 600 \
  --round both --title "Amino acid composition comparison" \
  --subtitle "Distribution across de novo and existing binders"
```

More R plotting scripts
```bash
# Density plots by category (e.g., metric distribution by round)
Rscript scripts/plotting_r/density.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/densities \
  --metric sequence_length \
  --category round \
  --format svg --width 1600 --height 1200 --res 300 \
  --round both

# Pairwise correlations with colored groups
Rscript scripts/plotting_r/correlation_plots.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/correlations \
  --x_column iptm \
  --y_column kd \
  --color_by design_category \
  --round both \
  --format svg --width 4000 --height 3600 --res 600

# Violin plots (e.g., KD across rounds)
Rscript scripts/plotting_r/violin_plots.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/violin \
  --y_column kd \
  --x_column round \
  --color_by round \
  --round both \
  --format svg --width 1600 --height 1200 --res 300 \
  --show_anova TRUE --binders_only FALSE

# Binding affinity ordered scatter with references
Rscript scripts/plotting_r/binding_affinity_plot.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/binding_affinity \
  --format svg --width 6000 --height 4000 --res 600

# 2x2 combined correlations vs KD
Rscript scripts/plotting_r/combined_metrics_plot.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/combined_metrics \
  --y_column kd \
  --color_by design_category \
  --round both \
  --format png --width 8000 --height 6000 --res 600 \
  --main_title "Correlation of Protein Metrics with Binding Affinity"

# Barplots for model types per round
Rscript scripts/plotting_r/barplots_model_types.R \
  --input ./data/processed/all_submissions_new.csv \
  --x_column round \
  --color_column RFdiffusion \
  --output ./plots/barplots_model_types \
  --format svg \
  --width 1600 --height 1400 --res 300 \
  --round both --title "Expressed" \
  --top_n 10 --subtitle "Number of expressed designs per round" \
  --sort alpha

# Interface property violins (multiple metrics)
Rscript scripts/plotting_r/violin_properties.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/interface_violins \
  --format svg --width 5000 --height 3800 --res 600 \
  --round both

# Radar plots for top binders
Rscript scripts/plotting_r/radar_plot.R \
  --input ./data/processed/all_submissions.csv \
  --output ./plots/radar_plots \
  --format png --width 3600 --height 3600 --res 600 \
  --round 2 \
  --binder_type both \
  --top_n 5 \
  --title "Interface metrics for the top binders" \
  --subtitle "Comparing De novo and Existing binders"
```

Binding property annotation pipeline
Two components exist: structure/interface annotation and ESM PLL scoring. Both are defined as Modal functions.

1) Interface/binding properties (Modal)
- Requirements are captured in `annotate_binding_properties.py` via a Modal `Image` builder. It installs utilities, sets up PyRosetta, and clones BindCraft.
- Customize local directories if needed:
  - `data/raw/structures/001_2024`, `002_2024`
  - `data/processed`
```bash
python annotate_binding_properties.py
```
This will batch process submissions, relax structures (if enabled), compute interface metrics, and produce CSVs under `data/processed`/Modal volume outputs.

2) ESM PLL scoring (Modal, GPU)
```bash
python annotate_esm_pll.py
```
This reads `data/processed/all_submissions.csv` (or the mounted path inside Modal), computes PLL scores, and writes results under the Modal volume.



