#!/usr/bin/env Rscript

# ------------------------------------------------------------
# SNP Gene-Track Plot (no haplotype blocks, no LD)
# ------------------------------------------------------------
# Usage:
#   Rscript snp_gene_track_plot.R <chromosome> <position> [WINDOW_bp] [GFF_path] [OUT_PREFIX]
# Example:
#   Rscript snp_gene_track_plot.R CM081021.1 47255812 100000 genomic.gff
#
# What it does:
#   • Reads a GFF/GFF3 file
#   • Extracts genes in a window around the SNP
#   • Draws gene boxes with strand arrows and labels
#   • Marks the SNP with a red dashed line
#   • Saves PDF/PNG and a TSV with genes in the region
# ------------------------------------------------------------

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(rtracklayer)
  library(scales)
  # ggrepel is optional; left out for portability
})

# -------------------------
# Parse args + defaults
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript snp_gene_track_plot.R <chromosome> <position> [WINDOW_bp] [GFF_path] [OUT_PREFIX]")
}

chrom   <- args[1]
pos     <- as.integer(args[2])
WINDOW  <- if (length(args) >= 3) as.integer(args[3]) else 100000L  # ±100 kb by default
GFF     <- if (length(args) >= 4) args[4] else "genomic.gff"
out_pref<- if (length(args) >= 5) args[5] else paste0(chrom, "_", pos, "_genetrack")

# -------------------------
# Region and helpers
# -------------------------
# Initial window (will be clipped to chromosome bounds below)
raw_start <- max(pos - WINDOW, 1L)
raw_end   <- pos + WINDOW

# Optional chromosome label mapping for prettier names in titles
chrom_map <- c(
  "CM081008.1" = "B1", "CM081009.1" = "B2", "CM081010.1" = "B3",
  "CM081011.1" = "B4", "CM081012.1" = "B5", "CM081013.1" = "B6",
  "CM081014.1" = "B7", "CM081015.1" = "B8", "CM081016.1" = "C1",
  "CM081017.1" = "C2", "CM081018.1" = "C3", "CM081019.1" = "C4",
  "CM081020.1" = "C5", "CM081021.1" = "C6", "CM081022.1" = "C7",
  "CM081023.1" = "C8", "CM081024.1" = "C9"
)
chrom_label <- if (chrom %in% names(chrom_map)) chrom_map[[chrom]] else chrom

# Chromosome sizes (bp); used to clip the window so it doesn't extend past ends
chrom_sizes <- c(
  "CM081008.1"=55420319, "CM081009.1"=71510277, "CM081010.1"=58930152,
  "CM081011.1"=58383510, "CM081012.1"=68377818, "CM081013.1"=64024298,
  "CM081014.1"=58421822, "CM081015.1"=69292429, "CM081016.1"=50873107,
  "CM081017.1"=62375533, "CM081018.1"=74943771, "CM081019.1"=64261724,
  "CM081020.1"=56455528, "CM081021.1"=47263926, "CM081022.1"=56217565,
  "CM081023.1"=49290977, "CM081024.1"=64281589
)

chr_len <- if (chrom %in% names(chrom_sizes)) as.integer(chrom_sizes[[chrom]]) else NA_integer_
if (!is.na(chr_len) && pos > chr_len) {
  stop(sprintf("SNP position (%s) exceeds chromosome length (%s)", format(pos, big.mark=","), format(chr_len, big.mark=",")))
}

if (!is.na(chr_len)) {
  region_start <- raw_start
  region_end   <- min(raw_end, chr_len)
  if (raw_start < 1L) region_start <- 1L
  if (raw_end > chr_len) message("Window clipped at chromosome end (", chrom_label, ": ", format(chr_len, big.mark=","), ")")
} else {
  # Fallback if size unknown
  region_start <- raw_start
  region_end   <- raw_end
}

region_id <- paste0(chrom_label, ":", format(pos, big.mark=","), " (±", WINDOW/1000, " kb)")

# Little helper: assign non-overlapping lanes to gene intervals
assign_lanes <- function(starts, ends, pad = 0L) {
  ord <- order(starts, ends)
  lanes_end <- integer(0)
  lane <- integer(length(starts))
  for (i in ord) {
    placed <- FALSE
    for (l in seq_along(lanes_end)) {
      if (starts[i] > (lanes_end[l] + pad)) {
        lane[i] <- l
        lanes_end[l] <- ends[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lanes_end <- c(lanes_end, ends[i])
      lane[i] <- length(lanes_end)
    }
  }
  lane
}

# Pick a reasonable display label per gene from available attributes
pick_gene_labels <- function(gr) {
  # potential attribute columns in order of preference
  prefs <- c("Name", "gene_name", "gene", "ID", "locus_tag", "Alias")
  have <- intersect(prefs, colnames(mcols(gr)))
  if (length(have) == 0) return(rep("gene", length(gr)))
  mat <- do.call(cbind, lapply(have, function(nm) as.character(mcols(gr)[[nm]])))
  # for each row, take first non-NA, non-empty
  apply(mat, 1, function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) "gene" else x[1]
  })
}

# -------------------------
# Read GFF and extract genes in window
# -------------------------
message("Reading GFF: ", GFF)
# Import and subset. If GFF is large, this still works; for bgz/indexed GFF it's efficient.
all_feats <- import(GFF)
if (!chrom %in% as.character(seqnames(all_feats))) {
  stop("Chromosome not found in GFF: ", chrom)
}

genes_chr <- all_feats[all_feats$type == "gene" & seqnames(all_feats) == chrom]
if (length(genes_chr) == 0) {
  warning("No 'gene' features found for chromosome ", chrom, ". Check your GFF.")
}

roi <- GRanges(seqnames = chrom, ranges = IRanges(start = region_start, end = region_end))
roi_genes <- subsetByOverlaps(genes_chr, roi)

# Prepare data.frame for plotting
if (length(roi_genes) > 0) {
  gene_df <- data.frame(
    start   = start(roi_genes),
    end     = end(roi_genes),
    width   = width(roi_genes),
    strand  = as.character(strand(roi_genes)),
    stringsAsFactors = FALSE
  )
  gene_df$label <- pick_gene_labels(roi_genes)
  # Lanes to avoid overlaps
  gene_df$lane <- assign_lanes(gene_df$start, gene_df$end, pad = 0L)
  # Strand-safe arrow endpoints at lane center
  gene_df$y <- gene_df$lane
  gene_df$x_start <- ifelse(gene_df$strand == "+", gene_df$start, gene_df$end)
  gene_df$x_end   <- ifelse(gene_df$strand == "+", gene_df$end, gene_df$start)
  nlanes <- max(gene_df$lane)
} else {
  gene_df <- data.frame(start=numeric(0), end=numeric(0), width=numeric(0),
                        strand=character(0), label=character(0), lane=integer(0),
                        y=numeric(0), x_start=numeric(0), x_end=numeric(0))
  nlanes <- 1
}

# -------------------------
# TSV summary of genes in region
# -------------------------
if (nrow(gene_df) > 0) {
  gene_df$overlaps_SNP <- pos >= gene_df$start & pos <= gene_df$end
  gene_df$dist_to_SNP  <- ifelse(gene_df$overlaps_SNP, 0L,
                                 pmin(abs(gene_df$start - pos), abs(gene_df$end - pos)))
  tsv <- gene_df[, c("label","start","end","strand","overlaps_SNP","dist_to_SNP")]
  setDT(tsv)
  setnames(tsv, "label", "gene")
  fwrite(tsv, paste0(out_pref, "_genes.tsv"), sep = "\t")
} else {
  fwrite(data.table(note = "No genes in window"), paste0(out_pref, "_genes.tsv"), sep = "\t")
}

# -------------------------
# Build plot
# -------------------------
p <- ggplot() +
  # Gene rectangles
  {if (nrow(gene_df) > 0) geom_rect(data = gene_df,
    aes(xmin = start, xmax = end,
        ymin = y - 0.3, ymax = y + 0.3,
        fill = strand), alpha = 0.6, color = NA)} +
  # Directional arrows along gene body
  {if (nrow(gene_df) > 0) geom_segment(data = gene_df,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed", ends = "last"))} +
  # Gene labels (at midpoint)
  {if (nrow(gene_df) > 0) geom_text(data = gene_df,
    aes(x = (start + end)/2, y = y + 0.45, label = label),
    size = 3, angle = 45, hjust = 0, vjust = 0.5, check_overlap = TRUE)} +
  # SNP vertical line
  geom_vline(xintercept = pos, linetype = "dashed", linewidth = 0.9, color = "red") +
  annotate("label", x = pos, y = nlanes + 0.9, label = paste0("SNP\n", pos),
           size = 3, label.size = 0.2, alpha = 0.9) +
  # Region frame
  coord_cartesian(xlim = c(region_start, region_end), ylim = c(0.3, max(1.5, nlanes + 1.2))) +
  scale_x_continuous(labels = label_number(si_suffix = TRUE)) +
  scale_fill_manual(values = c("+" = "#69b3a2", "-" = "#8da0cb"),
                    breaks = c("+","-"), labels = c("+ strand","- strand")) +
  labs(title = paste0("Gene track around ", region_id),
     x = "Genomic position (bp)", y = NULL, fill = "Strand") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# -------------------------
# Save outputs
# -------------------------
height <- max(3, nlanes * 0.8 + 2)

pdf_file <- paste0(out_pref, ".pdf")
png_file <- paste0(out_pref, ".png")

ggsave(pdf_file, plot = p, width = 16, height = height, units = "in")
ggsave(png_file, plot = p, width = 16, height = height, units = "in", dpi = 300)

message("✅ Saved: ", pdf_file)
message("✅ Saved: ", png_file)
message("✅ Saved: ", paste0(out_pref, "_genes.tsv"))
