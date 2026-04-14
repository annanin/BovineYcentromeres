# Bovine Sex Chromosome Centromere Analysis

**Author:** Anna Ní Nualláin
**Date:** April 2026

This repository contains the analysis pipeline and visualisation scripts for characterising the Y and X chromosome centromeres across six cattle breed assemblies. The Y centromere analysis is validated against the published Wagyu telomere-to-telomere assembly (Olagunju et al. 2024). The X centromere repeat analysis uses the five repeat sequences (XCTR1–5) characterised by Pineda et al. (2025).

├── cattle_Y_centromere_pipeline.sh     # Y chromosome centromere pipeline
├── 04_full_chromosome_distribution.sh  # Full Y and X chromosome monomer distribution
├── cattle_X_centromere_pipeline.sh     # X chromosome XCTR repeat pipeline 
├── 01_Y_data_preparation.R             # Load and prepare Y centromere data
├── 02_Y_figures.R                      # Generate all Y centromere figures
├── 03_X_figures.R                      # Generate all X centromere figures
└── README.md

---

## Dependencies

### Bash
- `blastn` / `makeblastdb` (BLAST+ v2.12.0 or later)
- `samtools` (v1.13 or later)
- `bedtools` (v2.30 or later)
- `trf` (Tandem Repeat Finder v4.09)
- `python3` (v3.8 or later; standard library only)
- `awk`, `bash`

### R
- `tidyverse` (v2.0 or later)
- `patchwork`
- `ggridges`

All bash scripts were developed and tested in a WSL2/Ubuntu environment on Windows.

---

## Input Files

Place the following files in `raw_data/` before running the Y pipeline:

| File | Contents |
|------|----------|
| `73bp_monomer.fasta` | 73 bp centromeric satellite monomer (Olagunju et al. 2024) |
| `Y_HOR.fasta` | Wagyu HOR consensus sequence |
| `Charolais_X-Wagyu_Y.fasta` | Wagyu Y chromosome (NC_082638.1) |
| `GCA_040286185.1_UOA_Wagyu_1_Tuli-Y.fasta` | Tuli Y chromosome (CM079978.1) |
| `Gyr_XY.fasta` | Gyr XY assembly |
| `Holstein_assembly.mixed-hap_mt.fasta` | Holstein XY assembly |
| `RedWagyu_plus-RubiaGallega-Y-unplaced.fasta` | Rubia Gallega Y chromosome |
| `Ayrshire_XY.fasta` | Ayrshire XY assembly |

Place the following in `raw_data/` before running the X pipeline:

| File | Contents |
|------|----------|
| `X_repeat_sequences.fasta` | XCTR1–5 sequences (Pineda et al. 2025) |

The same raw assembly FASTAs listed above are also used to extract X chromosomes in the X pipeline.

---

## Assemblies

Several assemblies are chimeric — the Y chromosome derives from a different individual than the X or autosomal sequence. All analyses refer to the Y or X chromosome breed as indicated.

| Y breed | X chromosome | Source assembly |
|---------|-------------|-----------------|
| Wagyu | Charolais | Charolais X + Wagyu Y (joined) |
| Tuli | Wagyu | GCA_040286185.1 (trio-phased) |
| Gyr | Gyr | Gyr XY assembly |
| Holstein | Holstein | Holstein mixed-haplotype assembly |
| Ayrshire | Ayrshire | Ayrshire XY assembly |
| Rubia Gallega | Red Wagyu | RedWagyu X + Rubia Gallega Y (joined) |

**Note on Ayrshire:** The Ayrshire assembly is in reverse complement orientation relative to the other five breeds. Centromere coordinates and monomer hit positions are corrected for this in both the bash and R scripts. The raw centromere position (77.3% along the Y) corrects to 22.7%, consistent with all other breeds (21.4–32.2%).

---

## Y Chromosome Pipeline

### Overview
cattle_Y_centromere_pipeline.sh
Step 1  Extract Y chromosomes from raw assemblies
Step 2  Extract candidate centromere regions (BLAST + bedtools merge)
Step 3  Monomer decomposition across the centromere
Step 4  TRF de novo structural assembly + HOR copy number by BLAST
04_full_chromosome_distribution.sh
Step 1  BLAST monomer against each full Y chromosome (>=60% identity)
Step 2  BLAST monomer against each X chromosome (>=60% identity)

### Running the pipeline

```bash
# Edit BASE_DIR in the script header to match your working directory, then:
bash cattle_Y_centromere_pipeline.sh
bash 04_full_chromosome_distribution.sh
```

Individual steps can be skipped by setting the corresponding `RUN_STEP_*` flag to `0` at the top of each script.

### Key methodological notes

**Centromere identification (Step 2):** The centromere is defined as the largest continuous cluster of merged BLAST hits from the 73 bp monomer and HOR consensus. One region per breed is selected. For Tuli, a secondary high-identity HOR block is visible outside the primary centromere boundaries at approximately 20–25 Mb on the full Y; without epigenetic data (CENP-A CUT&RUN or methylation) for this breed, it is not possible to determine whether this block is also a part of the active centromere.

**HOR copy number (Step 4):** Copy number is estimated by BLASTing the Wagyu HOR consensus at ≥90% identity and ≥80% query coverage against each full Y chromosome and counting non-overlapping hits within centromere coordinates. The 90% threshold was chosen empirically: the 95–98% identity range reported by Olagunju et al. (2024) refers to pairwise inter-copy identity within the array, not to identity between each copy and the consensus sequence. At ≥95%, only 228/364 Wagyu copies are recovered. At ≥90%, 380 copies are recovered (4% difference from published), which was considered satisfactory for cross-breed comparison. 

**TRF structural validation (Step 4):** TRF (v4.09) is run on each full Y assembly using the parameters of Olagunju et al. (2024). Two alternating repeat block types are identified de novo within the centromere, and inter-block spacer distances are measured without reference to published values. The Wagyu HOR start-to-start period recovered by TRF (3,740 bp) agrees closely with the published genomic repeat period (3,741 bp). De novo spacer distances (1,551–1,596 bp across all breeds) are consistently larger than the published satellite sequence spacer lengths (Spacer-1 = 1,478 bp; Spacer-2 = 1,460 bp), which likely reflects a difference in how TRF models the boundaries of diverged repeat blocks.

### Outputs

| Path | Contents |
|------|----------|
| `data/Ychromosomes/` | Extracted Y chromosome FASTAs per breed |
| `data/centromere_sequences/` | Extracted centromere FASTAs per breed |
| `analysis/monomer_decomposition/tsv/` | Per-breed annotated monomer hit tables |
| `analysis/monomer_decomposition/summary/` | Cross-breed monomer summary |
| `analysis/TRF/HOR_structural_results.tsv` | De novo TRF spacer distances per breed |
| `analysis/TRF/HOR_copies_BLAST.tsv` | BLAST-validated HOR copy numbers |
| `analysis/full_chromosome/tsv/` | Full Y and X monomer hit tables |
| `analysis/full_chromosome/summary/` | Centromere coordinates and X hit counts |

---

## X Chromosome Pipeline

### Overview
cattle_X_centromere_pipeline.sh
Step 1  Extract X chromosomes from raw assemblies
Step 2  Build BLAST databases
Step 3  BLAST XCTR1–5 at strict thresholds (>=90% id, >=80% cov)
Step 4  Copy number and centromere position summary
Step 5  Relaxed search for Gyr XCTR3 (absent at strict thresholds)
Step 6  Threshold sensitivity analysis for XCTR1 and XCTR2

### Running the pipeline

```bash
bash cattle_X_centromere_pipeline.sh
```

### Key methodological notes

**BLAST thresholds:** XCTR1 and XCTR2 have a substantial population of genuine copies with identity just below the 90% strict threshold. For these two repeats, relaxed threshold results (≥70% identity, ≥50% query coverage) are used in the final copy number analysis. XCTR3–5 are robust to threshold choice and are reported from strict threshold searches.

**Assembly orientation:** Ayrshire, Holstein, and Gyr X assemblies are in reverse complement orientation. Following coordinate correction, all six breeds localise centromeric repeats to 38–49 Mb from the p-arm telomere, consistent with Pineda et al. (2025).

**Gyr XCTR3:** XCTR3 (318 kb) is absent from Gyr at strict thresholds. A relaxed search recovers a partial copy (220 kb of 318 kb; 99.5% sequence identity) at approximately 42.4 Mb, suggesting a structural rearrangement or assembly truncation rather than sequence divergence.

### Outputs

| Path | Contents |
|------|----------|
| `data/*_X.fasta` | Extracted X chromosome FASTAs |
| `blast_results/*_XCTR_strict.txt` | Strict threshold BLAST hits per breed |
| `blast_results/*_XCTR_relaxed.txt` | Relaxed threshold BLAST hits per breed |
| `blast_results/Gyr_XCTR3_relaxed.txt` | Gyr XCTR3 relaxed search results |

---

## R Analysis Scripts

Run the scripts in order. All three assume the bash pipelines have been run first.

```r
# In R or RStudio:
source("01_Y_data_preparation.R")   # loads and cleans all Y data
source("02_Y_figures.R")            # generates Y centromere figures
source("03_X_figures.R")            # generates X centromere figures
```

**`01_Y_data_preparation.R`** loads all BLAST and TRF output tables, computes inter-monomer spacings, and assembles the final copy number estimates. It must be sourced before `02_Y_figures.R`, or the saved CSVs in `analysis/R_tables/` can be loaded manually.

**`02_Y_figures.R`** generates the following figures:

| File | Description |
|------|-------------|
| `Fig1_HOR_structure_schematic.pdf` | HOR unit schematic, full Y monomer density tracks, zoomed centromere HOR identity, copy number and array span |
| `Fig3_monomer_spacing.pdf` | Inter-monomer spacing histogram and category proportions |
| `Fig4_chromosomal_context.pdf` | Wagyu full Y monomer distribution and Y vs X specificity |
| `SuppFig1_monomer_position_tracks.pdf` | Per-breed monomer identity along the centromere |
| `SuppFig2_monomer_identity_distributions.pdf` | Monomer identity ridgeline plot |
| `SuppFig3_period_sensitivity.pdf` | Copy number sensitivity to assumed HOR period |
| `SuppFig4_structural_stability.pdf` | Structural Z-score summary across breeds |
| `SuppFig5_TRF_spacer_distances.pdf` | De novo TRF spacer distances vs published values |

**`03_X_figures.R`** generates the following figures:

| File | Description |
|------|-------------|
| `FigX_XCTR_full_chromosome.pdf` | XCTR1–5 hit tracks across the full X chromosome |
| `FigX_XCTR_centromere_zoom.pdf` | Zoomed centromeric region (30–55 Mb) |
| `FigX_XCTR_threshold_comparison.pdf` | Strict vs relaxed threshold comparison |

---

## References

Olagunju et al. (2024). Telomere-to-telomere assembly of the cattle sex chromosomes reveals extensive centromeric repeat variation. *[journal]*.

Pineda et al. (2025). Insights into the bovine X chromosome centromere structure and organisation. *[journal]*.

Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research*, 27(2), 573–580.

Camacho et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.
