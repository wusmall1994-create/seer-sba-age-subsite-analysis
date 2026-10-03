#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("SBA_FIGURE_DIR")
if (!nzchar(out_dir)) stop("Supply the figure output directory as the first argument or SBA_FIGURE_DIR.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(out_dir, "FigureS3_detailed_cohort_flow")

flow <- data.frame(
  step = c(
    "Master descriptive cohort",
    "Stage cohort",
    "Diagnosis window for survival analysis",
    "Primary survival cohort"
  ),
  detail = c(
    "14,158 adults with malignant adenocarcinoma-only SBA\nSEER 17, 2000-2023",
    "12,391 cases diagnosed in 2004-2023\nClassifiable stage4, including 1,049 unknown/unstaged\nExcluded: 1,767 diagnoses before 2004 or records without classifiable stage4",
    "11,592 cases diagnosed in 2004-2022\nExcluded: 799 diagnoses in 2023",
    "11,563 cases with valid survival follow-up\nExcluded: 29 DCO/autopsy cases without valid survival information\nIncludes 959 unknown/unstaged cases (8.29%)"
  ),
  stringsAsFactors = FALSE
)
write.csv(flow, paste0(stem, "_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

draw_flow <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 1), yscale = c(0, 1)))

  ys <- c(0.86, 0.62, 0.38, 0.14)
  fills <- c("#D9EAF7", "#E2F0D9", "#FFF2CC", "#FCE4D6")
  heights <- c(0.14, 0.18, 0.14, 0.18)

  for (i in seq_along(ys)) {
    grid::grid.roundrect(
      x = 0.50, y = ys[i], width = 0.76, height = heights[i], r = grid::unit(0.018, "snpc"),
      gp = grid::gpar(fill = fills[i], col = "#4D4D4D", lwd = 1.2)
    )
    grid::grid.text(
      flow$step[i], x = 0.50, y = ys[i] + heights[i] * 0.22,
      gp = grid::gpar(fontfamily = "Arial", fontsize = 10.5, fontface = "bold", col = "#1F1F1F")
    )
    grid::grid.text(
      flow$detail[i], x = 0.50, y = ys[i] - heights[i] * 0.09,
      gp = grid::gpar(fontfamily = "Arial", fontsize = 8.4, lineheight = 1.05, col = "#1F1F1F")
    )
    if (i < length(ys)) {
      grid::grid.lines(
        x = c(0.50, 0.50),
        y = c(ys[i] - heights[i] / 2 - 0.010, ys[i + 1] + heights[i + 1] / 2 + 0.018),
        arrow = grid::arrow(type = "closed", length = grid::unit(0.13, "inches")),
        gp = grid::gpar(col = "#4D4D4D", lwd = 1.2)
      )
    }
  }
  grid::popViewport()
}

width_in <- 7.2
height_in <- 8.0

grDevices::png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 600, bg = "white", type = "cairo")
draw_flow()
grDevices::dev.off()

grDevices::tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", bg = "white", type = "cairo")
draw_flow()
grDevices::dev.off()

grDevices::svg(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white", family = "Arial")
draw_flow()
grDevices::dev.off()

grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial", bg = "white")
draw_flow()
grDevices::dev.off()

message("Created detailed cohort flow exports at: ", stem)
