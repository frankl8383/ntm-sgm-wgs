#!/usr/bin/env Rscript

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(paste("Unexpected argument:", key))
    }
    if (i == length(args)) {
      stop(paste("Missing value for", key))
    }
    parsed[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  parsed
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "summary", "pdf", "png", "top-n", "width", "height")
missing <- required[!required %in% names(args)]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

suppressPackageStartupMessages(library(ggplot2))

species_long <- read_tsv(args[["input"]])
summary <- read_tsv(args[["summary"]])
top_n <- as.integer(args[["top-n"]])
plot_width <- as.numeric(args[["width"]])
plot_height <- as.numeric(args[["height"]])

if (nrow(species_long) == 0) {
  plot_data <- data.frame(
    sample_id = factor(summary$sample_id, levels = summary$sample_id),
    species_group = "No species assignment",
    fraction_total_reads = 1
  )
} else {
  species_long$fraction_total_reads <- as.numeric(species_long$fraction_total_reads)
  species_totals <- aggregate(fraction_total_reads ~ species, species_long, sum)
  species_totals <- species_totals[order(species_totals$fraction_total_reads, decreasing = TRUE), , drop = FALSE]
  top_species <- head(species_totals$species, top_n)

  plot_data <- species_long
  plot_data$species_group <- ifelse(plot_data$species %in% top_species, plot_data$species, "Other")
  plot_data <- aggregate(
    fraction_total_reads ~ sample_id + species_group,
    plot_data,
    sum
  )

  observed <- split(plot_data$fraction_total_reads, plot_data$sample_id)
  missing_rows <- list()
  for (sample_id in summary$sample_id) {
    total <- sum(observed[[sample_id]], na.rm = TRUE)
    remainder <- max(0, 1 - total)
    if (remainder > 1e-6) {
      missing_rows[[length(missing_rows) + 1]] <- data.frame(
        sample_id = sample_id,
        species_group = "Unassigned",
        fraction_total_reads = remainder
      )
    }
  }
  if (length(missing_rows) > 0) {
    plot_data <- rbind(plot_data, do.call(rbind, missing_rows))
  }
  plot_data$sample_id <- factor(plot_data$sample_id, levels = summary$sample_id)
}

if (exists("top_species")) {
  species_levels <- c(top_species[top_species %in% plot_data$species_group], "Other", "Unassigned")
  species_levels <- species_levels[species_levels %in% plot_data$species_group]
} else {
  species_levels <- unique(plot_data$species_group)
}
plot_data$species_group <- factor(plot_data$species_group, levels = species_levels)

p <- ggplot(plot_data, aes(x = sample_id, y = fraction_total_reads, fill = species_group)) +
  geom_col(width = 0.82, color = "grey20", linewidth = 0.08) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), expand = c(0, 0)) +
  labs(x = NULL, y = "Read fraction", fill = "Species") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.key.height = unit(0.45, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )

dir.create(dirname(args[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[["pdf"]], p, width = plot_width, height = plot_height, units = "in")
ggsave(args[["png"]], p, width = plot_width, height = plot_height, units = "in", dpi = 300)
