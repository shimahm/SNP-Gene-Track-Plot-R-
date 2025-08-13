# SNP Gene-Track Plot (R)

A small, dependency-light R script to visualize **genes around a focal SNP** on a given chromosome.  
It draws strand-aware gene boxes with arrows, marks the SNP with a dashed red line, and exports high-quality PDF/PNG plus a TSV summary of genes in the window.

> No LD heatmaps, no haplotype blocks — just a clean gene track around your SNP.

---

## Features

- **Input:** chromosome ID, SNP position (bp), optional window size and GFF path
- **Output:**
  - Publication-ready **PDF** and **PNG**
  - **TSV** listing genes in the window (whether they overlap the SNP and distance to SNP)
- Strand arrows and non-overlapping lanes for dense regions
- Pretty chromosome labels (e.g., `CM081021.1 → C6`)
- **Window automatically clipped** to real chromosome length (prevents plotting beyond contig ends)

---

## Requirements

- **R ≥ 4.1**
- Packages:
  - `data.table`
  - `ggplot2`
  - `GenomicRanges` (Bioconductor)
  - `rtracklayer` (Bioconductor)
  - `scales`

### Install packages

```r
install.packages(c("data.table", "ggplot2", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

---

## Usage

```bash
Rscript snp_gene_track_plot.R <chromosome> <position> [WINDOW_bp] [GFF_path] [OUT_PREFIX]
```

- `<chromosome>`: contig name in your GFF (e.g., CM081021.1)
- `<position>`: 1-based genomic coordinate of the SNP (e.g., 47255812)
- `[WINDOW_bp]` (optional): half-window in bp (default 100,000 → ±100 kb)
- `[GFF_path]` (optional): path to your GFF/GFF3 (default `genomic.gff`)
- `[OUT_PREFIX]` (optional): file prefix (default `<chromosome>_<position>_genetrack`)

### Examples

```bash
# ±100 kb around SNP on C6
Rscript snp_gene_track_plot.R CM081021.1 47255812 100000 genomic.gff

# ±250 kb, custom output prefix
Rscript snp_gene_track_plot.R CM081020.1 31630794 250000 genomic.gff C5_31630794
```

---

## Inputs

- GFF/GFF3 with `type == "gene"` rows for gene bodies.
- The script chooses a label from Name, gene_name, gene, ID, locus_tag, or Alias (first available).

---

## Outputs

For `OUT_PREFIX = C6_47255812_genetrack`:

- `C6_47255812_genetrack.pdf`
- `C6_47255812_genetrack.png`
- `C6_47255812_genetrack_genes.tsv` with columns:
  - `gene` (label)
  - `start`, `end`
  - `strand` (+ / -)
  - `overlaps_SNP` (TRUE/FALSE)
  - `dist_to_SNP` (0 if overlapping; otherwise nearest bp distance)

---

## Customization

### 1) Pretty chromosome labels

Defined near the top of the script:

```r
chrom_map <- c(
  "CM081008.1" = "B1", "CM081009.1" = "B2", "CM081010.1" = "B3",
  "CM081011.1" = "B4", "CM081012.1" = "B5", "CM081013.1" = "B6",
  "CM081014.1" = "B7", "CM081015.1" = "B8", "CM081016.1" = "C1",
  "CM081017.1" = "C2", "CM081018.1" = "C3", "CM081019.1" = "C4",
  "CM081020.1" = "C5", "CM081021.1" = "C6", "CM081022.1" = "C7",
  "CM081023.1" = "C8", "CM081024.1" = "C9"
)
```

### 2) Chromosome sizes (window clipping)

Provided in the script (bp):

```r
chrom_sizes <- c(
  "CM081008.1"=55420319, "CM081009.1"=71510277, "CM081010.1"=58930152,
  "CM081011.1"=58383510, "CM081012.1"=68377818, "CM081013.1"=64024298,
  "CM081014.1"=58421822, "CM081015.1"=69292429, "CM081016.1"=50873107,
  "CM081017.1"=62375533, "CM081018.1"=74943771, "CM081019.1"=64261724,
  "CM081020.1"=56455528, "CM081021.1"=47263926, "CM081022.1"=56217565,
  "CM081023.1"=49290977, "CM081024.1"=64281589
)
```
The plotting window is clipped to [1, chromosome_length]. If a SNP is beyond the chromosome end, the script stops with a clear error.

### 3) Figure size & legend

- Width is set to 16 in; height auto-scales with number of gene lanes.
- Legend is on the right.
- To change width, edit the two `ggsave(..., width = 16, height = height, ...)` lines.

---

## Repo Structure (suggested)

```
.
├── snp_gene_track_plot.R
├── data/
│   └── genomic.gff               # your annotation (not tracked if large)
├── examples/
│   └── example_plot.png          # add a thumbnail here
├── README.md
└── LICENSE
```

---

## Troubleshooting

- **“Chromosome not found in GFF”**  
  Ensure the chromosome argument matches seqnames in your GFF (e.g., `CM081021.1`), not the pretty label (C6).

- **“No 'gene' features found …”**  
  Your GFF may use different feature types. Adjust the filter: `all_feats$type == "gene"`.

- **Empty region**  
  If no genes fall in the window, the plot still renders with the SNP line; the TSV will contain a note.

---

## Reproducibility

You can embed your sessionInfo() to the README of a specific run:

```r
sessionInfo()
```

---

## Citation

If this script helps your research, please cite this repository and the R packages used:  
`data.table`, `ggplot2`, `GenomicRanges`, `rtracklayer`, `scales`.

---

## License

MIT

---


