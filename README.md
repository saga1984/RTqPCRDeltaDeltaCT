# RTqPCRDeltaDeltaCT

Utilities to process RT-qPCR data (Bio-Rad CFX96 exports) and produce clean, publication-ready results using the 2^-ΔΔCt (Livak) method. The workflow is split into two R scripts:

1. **`RTqPCR_data_analysis`** — cleans raw Cq data, trims replicates, computes ΔΔCt / Fold Change, and exports per-comparison tables + plots (including a **Cq boxplot to assess dispersion**).
2. **`RTqPCR_data_analysis_unir_datos`** — aggregates those per-comparison tables across treatments and reference-gene sets, builds combined bar/box plots with error bars and optional stats (ANOVA/Tukey), and exports summary files.

---

## Features

* Reads Bio-Rad CFX96 CSVs (`Target`, `Content`, `Sample`, `Cq`)
* Flags NTCs and missing Cq; trims to the closest *n* replicates per group
* Supports **single or multiple reference genes** (averaged Cq)
* Computes ΔCt, ΔΔCt (calibrator vs treatment), and **Fold Change = 2^-ΔΔCt**
* Exports:

  * Cleaned/filtered table per run
  * Per-comparison summary (`Average`, `SD`)
  * Aggregated tables across runs/refs
* **Plots**:

  * **Cq boxplot** (dispersion/quality check)
  * Bar plots (+ SE) for Fold Change and log2(Fold Change)
  * Boxplots (with jitter) for Fold Change and log2(Fold Change)
* Optional **ANOVA + Tukey HSD** and Excel export of stats
* Saves figures as TIFF/JPEG and optionally PDF/PNG/SVG

---

## Requirements

R ≥ 4.1 with packages:

```
psych, dplyr, ggplot2, tidyr, ggrepel, ggpubr, openxlsx, svglite
```

Install once:

```r
install.packages(c("psych","dplyr","ggplot2","tidyr","ggrepel","ggpubr","openxlsx","svglite"))
```

---

## Input data format

CSV exported from CFX96 with **exact** columns:

* `Target` — gene name
* `Content` — e.g., `NTC`, `Unknown`
* `Sample` — condition label (e.g., `0 Hours`, `15 Minutes`, `Treatment`)
* `Cq` — numeric Cq

> The script orders rows by `Target`, `Content`, `Sample`. NTCs with `Cq >= 31` are flagged as suspicious; rows with `NA(Cq)` are dropped from analysis.

---

## Script 1 — `RTqPCR_data_analysis`

### What it does

* Shows available **Targets** and **Samples**, then **interactively asks** for:

  * gene of interest (`gen_interes`)
  * reference gene(s) (comma or space separated)
  * calibrator condition (control)
  * treatment condition
* For multiple reference genes, averages their Cq before ΔCt.
* Trims to the closest `numero_replicas` (default 3) per `Target`×`Sample`.
* Computes ΔCt, ΔΔCt (vs calibrator), **Fold Change = 2^-ΔΔCt**.
* **Plots:**

  * **Cq boxplot** per `Target` (**dispersion/quality**)
  * Bar plot of Fold Change with SD error bars

### Example

```r
source("RTqPCR_data_analysis.R")

RTqPCR_data_analysis(
  ruta_rtqpcr = "2025/articulo_Rosy/target_genes.csv",
  numero_replicas = 3,
  resolucion = 600,
  formatos = "jpeg"   # or c("tiff","jpeg")
)
```

### Outputs (same folder as input)

* `Boxplot_<file>.{ext}` — **Cq boxplot (dispersion)**
* `Tabla_filtrada_<file>_<tratamiento>_<GOI>_<REFS>.csv` — filtered rows used
* `Tabla_Fold_Change_Data_<file>_<tratamiento>_<GOI>_<REFS>.csv` — per-replicate ΔCt/ΔΔCt/Fold Change
* `Fold_Change_Data_<file>_<tratamiento>_<GOI>_<REFS>.csv` — per group means/SD
* `fold_Livak_<file>_<tratamiento>_<GOI>_<REFS>.{ext}` — Fold Change bar plot

> **Note on naming:** with multiple reference genes, names are joined by `_` (e.g., `Ubiq_Rps1`).
> **Note on columns:** the script stores the Fold Change (2^-ΔΔCt) in a column named `log2FC`. If you truly need log2 scale, compute `log2(Fold Change)` downstream.

> [Boxplot](man/figures/boxplot_all.jpeg)

> [Fold Change Treatment vs Control] (fold_Livak_articulo_endogenos_Quantification_Cq_Results_1_18Jul2024_ABA_ABI5_PsbA_ATP6.jpeg)

---

## Script 2 — `RTqPCR_data_analysis_unir_datos`

### What it does

* Scans a directory for outputs from Script 1 and **binds** them:

  * Aggregated “bar” tables from `Fold_Change_Data_*`
  * Aggregated “box” tables from `Tabla_Fold_Change_Data_*`
* Creates combined **bar plots** (mean ± SE) and **boxplots** (with optional jitter)
* Allows **facet by reference-gene set** and adds significance brackets
* Optionally runs **ANOVA** (per reference set) and **Tukey HSD**, and writes an Excel workbook

### Minimal example

```r
source("RTqPCR_data_analysis_unir_datos.R")

RTqPCR_data_analysis_unir_datos(
  gen_interes   = "ERF03",
  condicion     = "Dehydration",
  gen_referencia= c("Ubiq","Rps1"),   # or "Ubiq"
  ruta_tablas   = "2025/articulo_Rosy/",
  formatos      = c("tiff","jpeg"),
  fig_ancho     = 5000,
  fig_alto      = 4000
)
```

### Additional combined plots & stats

The script includes examples to:

* Bind multiple `big_df_*` objects and facet by `Reference_Gene`
* Add significance brackets and **`ggpubr::stat_compare_means()`**
* Run **ANOVA/Tukey** per dataset and export:

  * `Resultados_ANOVA_Tukey.xlsx`
  * `stats_<GEN>.xlsx`

### Key outputs

* Aggregated CSVs:

  * `Fold_Change_agrupado_<condicion>_<GOI>_<REFS>.csv`
  * `Tabla_Fold_Change_agrupado_<condicion>_<GOI>_<REFS>.csv`
* Figures (saved to `ruta_tablas`):

  * `FC_barplot_agrupado_*.{tiff,jpeg}` and `log2FC_barplot_agrupado_*`
  * `FC_boxplot_agrupado_*` and `log2FC_boxplot_agrupado_*`
  * Optional `*.pdf`, `*.png`, `*.svg` via **svglite**
* **Global R objects** (for reuse):
  `big_df_bar_*`, `big_df_box_*`, `FC_barplot_agrupado_*`, `log2FC_barplot_agrupado_*`, etc., are assigned into `.GlobalEnv`.

> [FoldChange Grouped Plot](man/figures/FC_BoxPlot_Agrupados_Estadisticos_ERF03.jpeg)

---

## Typical workflow

1. **Prepare CSV** from CFX96 with the required columns.
2. **Run Script 1** once per gene/treatment/reference-set to generate per-comparison tables & plots.
3. **Run Script 2** to aggregate all those outputs for a given `gen_interes` and `condicion`, and to build combined **bar/box plots** (plus stats, if desired).
4. Use the exported **Cq boxplot** to assess dispersion and replicate consistency.
5. Use aggregated figures/tables for the manuscript or report.

---

## Plot notes

* **Cq boxplot** (Script 1): quick check of dispersion/outliers per Target.
* **Bar plots**: show mean Fold Change or log2(Fold Change) with **SE** error bars across `Treatment`.
* **Boxplots**: visualize **distribution/dispersion** of values (optionally overlaid jitter).
* Colors are fixed by sample (`Control`, `Treatment`) and can be edited in the code (`colores_fijos`).

---

## Tips & troubleshooting

* **Multiple reference genes**: their Cq values are averaged before ΔCt. Ensure names in the CSV exactly match your input.
* **Replicate trimming**: set `numero_replicas` to the number you expect per group (default 3).
* **Factor order**: treatments are ordered as `0 Hours → 15 Minutes → 30 Minutes → 1 Hour → 3 Hours → 24 Hours`. Edit the factor levels if yours differ.
* **Large y-ranges**: adjust `ylim` or compute per-facet maxima to keep labels visible.
* **Stats**: `ggpubr` must be installed; ANOVA/Tukey require sufficient replicates.
* **Windows paths**: use forward slashes (`"2025/articulo_Rosy/"`) or wrap backslashes.

---

## Citation

If you use this workflow, please cite the original method:

* **Livak & Schmittgen (2001)**. Analysis of Relative Gene Expression Data Using Real-Time Quantitative PCR and the 2^-ΔΔCt Method. *Methods* 25(4):402–408.

---

## License

Specify your license (e.g., MIT) in `LICENSE`.

---

## Acknowledgments

Developed for CFX96 RT-qPCR analyses with emphasis on transparent QC (Cq dispersion), robust ΔΔCt computation, and publication-grade outputs.

## 👩‍💻 Author
Developed with ❤️ in R by [LinkedIn](https://www.linkedin.com/in/santiagovalentingalvangordillo) | [ORCID](https://orcid.org/0000-0001-6609-5661)   
