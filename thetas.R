# thetas.R (refactored) - simplified and de-duplicated
# Set working directory if needed:
# setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sfs/thetas")

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)

# Column names used in ANGSD pestPG output windows
colnames_vec <- c("metrics", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima",
                  "fuf", "fud", "fayh", "zeng", "nSites")

# Named vector mapping sample codes -> filenames
sample_files <- c(
  TYD = "TYD.thetas.windows.gz.pestPG",
  TYE = "TYE.thetas.windows.gz.pestPG",
  TYA = "TYA.thetas.windows.gz.pestPG",
  HAY = "HAY1.thetas.windows.gz.pestPG",
  CAM = "CAM.thetas.windows.gz.pestPG",
  CWA = "CWA.thetas.windows.gz.pestPG",
  CWB = "CWB.thetas.windows.gz.pestPG",
  CWC = "CWC.thetas.windows.gz.pestPG",
  WEA = "WEA.thetas.windows.gz.pestPG",
  WEF = "WEF.thetas.windows.gz.pestPG",
  AGA = "AGA.thetas.windows.gz.pestPG",
  AGB = "AGB.thetas.windows.gz.pestPG",
  AVB = "AVB.thetas.windows.gz.pestPG",
  AVD = "AVD.thetas.windows.gz.pestPG",
  OWE = "OWE.thetas.windows.gz.pestPG",
  NEV = "NEV2.thetas.windows.gz.pestPG",  # original file was NEV2 -> refer to it as NEV in code
  YHA = "YHA.thetas.windows.gz.pestPG",
  YHE = "YHE.thetas.windows.gz.pestPG"
)

# Helper to read a file and standardise column names and compute Pi
read_theta <- function(filename) {
  df <- read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df) <- colnames_vec
  df$Pi <- df$tP / df$nSites
  df
}

# Read all samples into a named list
thetas <- map(sample_files, read_theta)

# ---- Helper plotting functions ----

# Create a list of line plots for a given chromosome region and metric (Tajima or Pi)
make_chr_line_plots <- function(samples, chr = "NC_042981.1", metric = "Tajima", xlim = c(20000000, 30000000)) {
  plots <- map(samples, function(s) {
    dat <- thetas[[s]] %>% filter(Chr == chr)
    ggplot(dat) +
      geom_line(aes_string(x = "WinCenter", y = metric)) +
      scale_x_continuous(limits = xlim) +
      theme_bw() +
      labs(title = s)
  })
  names(plots) <- samples
  plots
}

# Create a density plot comparing multiple samples for a given metric
make_density_group_plot <- function(samples, metric = "Tajima", colors = NULL,
                                    xlim = NULL, ylim = NULL, legend = FALSE) {
  # build matrix/dataframe of metric columns
  metric_cols <- map_dfc(samples, ~ thetas[[.x]][[metric]])
  names(metric_cols) <- samples
  # remove rows containing any zero (matching original behaviour)
  metric_cols <- metric_cols[rowSums(metric_cols == 0) == 0, , drop = FALSE]
  long <- pivot_longer(metric_cols, cols = everything(), names_to = "sample", values_to = "value")
  p <- ggplot(long, aes(x = value, color = sample)) +
    geom_density(alpha = 0.5) +
    labs(x = ifelse(metric == "Pi", "Pi", "Tajima's D"), y = "Density") +
    theme_bw() +
    theme(legend.position = ifelse(legend, "right", "none"))
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)
  if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim)
  if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim)
  p
}

# ---- Reproduce the outputs from your original script ----

# 1) Chromosome (NC_042981.1) line plots (TYA, TYD, TYE, WEA, WEF) for Tajima and Pi
chr_samples <- c("TYA", "TYD", "TYE", "WEA", "WEF")
taj_plots <- make_chr_line_plots(chr_samples, metric = "Tajima", xlim = c(20000000, 30000000))
pi_plots  <- make_chr_line_plots(chr_samples, metric = "Pi",     xlim = c(20000000, 30000000))

# combine in the same order (as original): TYA, TYD, TYE, WEA, WEF (one column)
chr25_ne_taj <- plot_grid(plotlist = taj_plots, ncol = 1)
ggsave("chr25_ne_taj.png", chr25_ne_taj, width = 6, height = 10)

chr25_ne_Pi <- plot_grid(plotlist = pi_plots, ncol = 1)
ggsave("chr25_ne_Pi.png", chr25_ne_Pi, width = 6, height = 10)

# 2) Density comparisons for various groups (with same axis limits/colors as original)
taj_xlim <- c(-2, 3)
taj_ylim <- c(0, 0.7)
pi_xlim  <- c(0, 0.5)
pi_ylim  <- c(0, 5)

# Tyne group
Tyne_Taj <- make_density_group_plot(c("TYD", "TYE", "TYA"), metric = "Tajima",
                                    colors = c("TYD" = "green", "TYE" = "orange", "TYA" = "deeppink"),
                                    xlim = taj_xlim, ylim = taj_ylim)
Tyne_Pi  <- make_density_group_plot(c("TYD", "TYE", "TYA"), metric = "Pi",
                                    colors = c("TYD" = "green", "TYE" = "orange", "TYA" = "deeppink"),
                                    xlim = pi_xlim, ylim = pi_ylim)

# Kernow
Kernow_Taj <- make_density_group_plot(c("HAY", "CAM"), metric = "Tajima",
                                      colors = c("HAY" = "green", "CAM" = "orange"),
                                      xlim = taj_xlim, ylim = taj_ylim)
Kernow_Pi  <- make_density_group_plot(c("HAY", "CAM"), metric = "Pi",
                                      colors = c("HAY" = "green", "CAM" = "orange"),
                                      xlim = pi_xlim, ylim = pi_ylim)

# Conwy
Conwy_Taj <- make_density_group_plot(c("CWA", "CWB", "CWC"), metric = "Tajima",
                                     colors = c("CWA" = "orange", "CWB" = "deeppink", "CWC" = "green"),
                                     xlim = taj_xlim, ylim = taj_ylim)
Conwy_Pi  <- make_density_group_plot(c("CWA", "CWB", "CWC"), metric = "Pi",
                                     colors = c("CWA" = "orange", "CWB" = "deeppink", "CWC" = "green"),
                                     xlim = pi_xlim, ylim = pi_ylim)

# Ireland
Ireland_Taj <- make_density_group_plot(c("AVB", "AVD", "OWE"), metric = "Tajima",
                                       colors = c("AVB" = "orange", "AVD" = "deeppink", "OWE" = "green"),
                                       xlim = taj_xlim, ylim = taj_ylim)
Ireland_Pi  <- make_density_group_plot(c("AVB", "AVD", "OWE"), metric = "Pi",
                                       colors = c("AVB" = "orange", "AVD" = "deeppink", "OWE" = "green"),
                                       xlim = pi_xlim, ylim = pi_ylim)

# Anglesey
Anglesey_Taj <- make_density_group_plot(c("AGA", "AGB", "NEV"), metric = "Tajima",
                                        colors = c("AGA" = "orange", "AGB" = "deeppink", "NEV" = "green"),
                                        xlim = taj_xlim, ylim = taj_ylim)
Anglesey_Pi  <- make_density_group_plot(c("AGA", "AGB", "NEV"), metric = "Pi",
                                        colors = c("AGA" = "orange", "AGB" = "deeppink", "NEV" = "green"),
                                        xlim = pi_xlim, ylim = pi_ylim)

# Wear
Wear_Taj <- make_density_group_plot(c("WEA", "WEF"), metric = "Tajima",
                                    colors = c("WEA" = "orange", "WEF" = "green"),
                                    xlim = taj_xlim, ylim = taj_ylim)
Wear_Pi  <- make_density_group_plot(c("WEA", "WEF"), metric = "Pi",
                                    colors = c("WEA" = "orange", "WEF" = "green"),
                                    xlim = pi_xlim, ylim = pi_ylim)

# Ystwyth (note NEV included)
Ystwyth_Taj <- make_density_group_plot(c("NEV", "YHA", "YHE"), metric = "Tajima",
                                       colors = c("NEV" = "blue", "YHA" = "green", "YHE" = "orange"),
                                       xlim = taj_xlim, ylim = taj_ylim)
Ystwyth_Pi  <- make_density_group_plot(c("NEV", "YHA", "YHE"), metric = "Pi",
                                       colors = c("NEV" = "blue", "YHA" = "green", "YHE" = "orange"),
                                       xlim = pi_xlim, ylim = pi_ylim)

# Combine all population Tajima/Pi plots in similar layout to original
combined_nucs <- plot_grid(
  Anglesey_Taj, Anglesey_Pi,
  Ireland_Taj, Ireland_Pi,
  Kernow_Taj, Kernow_Pi,
  Tyne_Taj, Tyne_Pi,
  Wear_Taj, Wear_Pi,
  Conwy_Taj, Conwy_Pi,
  Ystwyth_Taj, Ystwyth_Pi,
  ncol = 2, rel_widths = c(3, 4)
)
ggsave("combined_pops_Taj_Pi.png", combined_nucs, width = 8, height = 14)

# Example quick single-plot view (as in your interactive checks)
# ggplot(thetas$TYD) + geom_density(aes(x = Pi)) + theme_bw()

# End of script
