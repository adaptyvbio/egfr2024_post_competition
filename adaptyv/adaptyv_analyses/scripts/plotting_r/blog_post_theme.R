library(ggplot2)
library(showtext)
library(ggtext)

font_add_google("Roboto", "Roboto")
showtext_auto()

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(linewidth = 0.5, color = "#E1E8ED"),
      panel.grid.minor = element_line(linewidth = 0.25, color = "#E1E8ED"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box.background = element_rect(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )
}

binding_colors <- c(
  "None" = "#E5E7EB",
  "Unknown" = "#F9F5C9", 
  "Weak" = "#C4B2FB",
  "Medium" = "#A7C1FB",
  "Strong" = "#2A4DD0"
)

expression_colors <- c(
  "None" = '#E5E7EB',
  "Low" = "#9EA2AF",
  "Medium" = "#E9C435", 
  "High" = "#69B7A7"
)

selection_status_colors <- c(
  "Top 100" = "#8B90DD",
  "Adaptyv selection" = "#8CD2F4",
  "Not selected" = "#3E6175"
)

design_category_colors <- c(
  "De novo" = "#56A6D4",
  "Optimized binder" = "#8B90DD",
  "Diversified binder" = "#8CD2F4",
  "Hallucination" = "#3E6175"
)

metric_colors <- c("ESM2 PLL" = "#8CD2F4", "ipTM" = "#8B90DD", "iPAE" = "#3E6175")

adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, selection_status_colors, design_category_colors, metric_colors)))
