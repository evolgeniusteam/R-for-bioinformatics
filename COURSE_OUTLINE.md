---
output:
  word_document: default
  html_document: default
---
# AI-powered bioinformatic analysis using R

**Instructor:** Wei-Hua Chen (HUST, China)
**License:** CC BY-NC 4.0
**Course Repository:** https://github.com/evolgeniusteam/R-for-bioinformatics

---

## Course Overview

16 lectures, divided into two parts:
- **Part 1 (Talks 1–6):** R Programming Foundations
- **Part 2 (Talks 7–16):** Data Analysis & Multi-omics Applications

---

## Part 1: R Programming Foundations

| Talk | Title | Source (merged from) |
|------|-------|---------------------|
| 01 | Course Overview & AI-assisted Setup | talk00 + talk01 + talk00_github |
| 02 | R Language Basics I | talk02 |
| 03 | R Language Basics II | talk03 |
| 04 | Data Wrangling with tidyverse | talk04 + talk05 |
| 05 | Data Visualization I — Base R & ggplot2 Fundamentals | talk06 (part) |
| 06 | Data Visualization II — Advanced ggplot2, Themes & Saving | talk06 (part) + talk09 (part) |
| 07 | Data Visualization III — Bioinformatics Plots | talk09 (part) + new content |
| 08 | Strings, Iteration & Functional Programming | talk07 + talk08 |

### Talk 01 — Course Overview & AI-assisted Setup
- **Course Overview**: about this course, 16-talk structure, target audience, learning approach
- **Why R?**: brief history (Ross Ihaka & Robert Gentleman, 1993), why R succeeded, R's popularity, R for bioinformatics
- AI-assisted programming philosophy with CodeBuddy
- **Step 1**: CodeBuddy installation & WeChat sign-in
- **Step 2**: Navigate CodeBuddy: interface, Agent modes (Ask / Craft / Plan)
- **Step 3**: Using prompts with Agent mode
- **Step 4**: Install R + RStudio via AI prompts
- **Step 5**: Install CodeBuddy extensions (R, R Debugger, Quarto)
- **Step 6–7**: Git clone course repo, open project in CodeBuddy + RStudio
- **Step 8**: Install tinytex & R packages (rmarkdown, knitr, bookdown, kableExtra)
- **Step 9**: Install tidyverse — the most important package
- **Step 10**: Render QMD to PDF — two verification exercises:
  - Verification 1: Beamer PDF (user name, R version, working directory, package list)
  - Verification 2: Homework-style PDF (same info, article format)
  - Manual rendering in RStudio (Render button)
- Verification checklist & prompt cheat sheet

### Talk 02 — R Language Basics I
- Data types: numeric, integer, character, logical, complex
- Vectors, matrices, arrays — creation & properties
- Vectorized operations (the key R paradigm)
- Basic arithmetic & logical operators
- Subsetting: `[ ]`, `[[ ]]`, `$`
- Recycling rules & common pitfalls
- Missing values (`NA`, `NaN`, `NULL`)
- Basic built-in functions (`length`, `sum`, `mean`, `summary`)

### Talk 03 — R Language Basics II
- Factors: creation, levels, ordering (`forcats`)
- Lists: creation, naming, subsetting
- Data frames vs tibbles — key differences
- Data import/export (`readr`: `read_csv`, `read_tsv`, `write_csv`)
- Introduction to tidy data concept
- Basic data frame operations (inspection, subsetting, sorting)
- Handling dates with `lubridate` (brief intro)

### Talk 04 — Data Wrangling with tidyverse
- The pipe operator (`%>%` / `|>`)
- `dplyr` core verbs:
  - `filter()` — subset rows by condition
  - `select()` — choose columns; helpers (`starts_with`, `contains`)
  - `mutate()` — create/modify columns; `across()` for column-wise ops
  - `arrange()` — sort rows
  - `summarise()` + `group_by()` — grouped summaries
- `tidyr` core functions:
  - `pivot_longer()` / `pivot_wider()` — reshape data
  - `separate()` / `unite()` — split/combine columns
- Joins: `left_join`, `inner_join`, `full_join`, `anti_join`
- `rowwise()` operations & nesting/unnesting

### Talk 05 — Data Visualization I — Base R & ggplot2 Fundamentals
- Why visualization matters in bioinformatics
- R base graphics:
  - `plot()`: scatter plots, parameters (type, pch, col, xlim, ylim)
  - `par()`: graphics parameters, `mfrow`/`mfcol` multi-panel
  - High-level vs. low-level plotting functions
  - Graphics devices: `pdf()`, `png()`, `dev.off()`
- Introduction to ggplot2:
  - Grammar of graphics philosophy
  - `ggplot()`, `aes()`, geom layers
  - Layers can use their own data
- Four core components: Layers, Scales, Coordinates, Facets
- Layers (geoms):
  - `geom_point()`, `geom_line()`, `geom_smooth()`
  - `geom_boxplot()`, `geom_violin()`
  - `geom_bar()` / `geom_col()`, `geom_histogram()`, `geom_density()`
- Scales: `scale_color_brewer()`, `scale_color_manual()`, `scale_color_gradient()`, `ggsci`
  - `colour` vs. `fill` distinction
- Coordinate systems: `coord_cartesian()`, `coord_flip()`, `coord_fixed()`, `coord_trans()`, `coord_polar()`, `coord_map()`
- Facets: `facet_grid()` / `facet_wrap()`

### Talk 06 — Data Visualization II — Advanced ggplot2, Themes & Saving
- Multi-panel plots:
  - `cowplot::plot_grid()` — grid layout
  - `ggdraw()` + `draw_plot()` — free-form layout
  - `gridExtra::grid.arrange()` — layout matrix
- Equations in plots:
  - `expression()` for static equations
  - `bquote()` / `substitute()` for variable substitution
  - Greek letters with `parse = TRUE`
  - Full example: regression equation on scatter plot
- Themes & legends:
  - Built-in themes: `theme_bw`, `theme_minimal`, `theme_classic`, etc.
  - Customizing with `theme()`: legend position, text size, grid lines
  - `labs()`: titles, axis labels, legend titles, tags
- Saving plots:
  - `ggsave()`: PDF, PNG, SVG with dpi control
  - Interactive plots: `plotly::ggplotly()`

### Talk 07 — Data Visualization III — Bioinformatics Plots with Real Data
- Common bioinformatics plot types overview
- Stacked bar plots for microbiome composition:
  - `geom_bar(stat = "identity", position = "stack")`
  - `reorder()` to control genus order
  - `position = "dodge"` for side-by-side comparison
  - Adding labels on stacked segments
  - Real data: Qin et al. gut microbiome genus data
- PCA for sample clustering:
  - `prcomp()` for dimensionality reduction
  - PCA scatter plot with variance explained
  - Interpretation: clusters, batch effects, outliers
- Box plots & violin plots for group comparisons:
  - Alpha diversity visualization (Shannon index, Richness)
  - Violin + box plot combo for distribution comparison
- Heatmaps with `geom_tile()`:
  - Log-transformation for abundance data
  - Color gradients for matrix visualization
  - Pattern interpretation (clusters, co-occurrence)
- Figure interpretation and storytelling in bioinformatics

### Talk 08 — Strings, Iteration & Functional Programming
- `stringr`: pattern matching, extraction, replacement
- Regular expressions: character classes, quantifiers, anchors, groups
- `purrr`:
  - `map()`, `map_dbl()`, `map_chr()`, `map_df()` — iterate & coerce output
  - `map2()`, `pmap()` — multi-input mapping
  - List-columns: storing complex objects in data frames
- Apply family: `lapply`, `sapply`, `vapply`, `tapply`
- Writing custom functions: arguments, return values, default values
- Defensive programming: `stopifnot()`, input validation
- (Brief) Parallel computing: `furrr` as `purrr` + `future`

---

## Part 2: Data Analysis & Multi-omics Applications

| Talk | Title |
|------|-------|
| 07 | Statistical Analysis & Machine Learning Basics |
| 08 | Regression & Survival Analysis |
| 09 | Transcriptome Analysis |
| 10 | Microbiome Analysis |
| 11 | Genome Analysis |
| 12 | Single-cell RNA-seq Analysis |
| 13 | Spatial Transcriptomics |
| 14 | Network Analysis |
| 15 | Multi-omics Integration & Case Studies |
| 16 | Other Omics & Course Wrap-up |

### Talk 07 — Statistical Analysis & Machine Learning Basics
- Descriptive statistics & summary tables
- Hypothesis testing: t-test, Wilcoxon, ANOVA, Kruskal-Wallis, χ²
- Multiple testing correction (Bonferroni, FDR/BH)
- Correlation analysis (Pearson, Spearman)
- Clustering: k-means, hierarchical clustering, dendrograms
- Dimensionality reduction: PCA, t-SNE, UMAP
- Classification basics: Random Forest
- Cross-validation & model evaluation

### Talk 08 — Regression & Survival Analysis
- Linear regression: model fitting, diagnostics, interpretation
- Logistic regression for binary outcomes
- Survival analysis basics: censoring, survival curves
- Kaplan-Meier estimator & log-rank test
- Cox proportional hazards model
- Clinical data preprocessing (missing values, encoding)
- Visualization: survival curves, forest plots, nomograms

### Talk 09 — Transcriptome Analysis
- RNA-seq experimental design & data structure
- Count matrix preprocessing
- Differential expression analysis: DESeq2 / edgeR / limma-voom
- Visualization: volcano plot, MA plot, heatmap, PCA plot
- Functional enrichment: GO (Gene Ontology), KEGG pathways
- GSEA (Gene Set Enrichment Analysis)
- Expression pattern clustering & visualization

### Talk 10 — Microbiome Analysis
- Microbiome data types: 16S rRNA amplicon, shotgun metagenomics
- OTU/ASV table, taxonomy table, metadata
- Alpha diversity: richness (Observed, Chao1), evenness (Shannon, Simpson)
- Beta diversity: Bray-Curtis, UniFrac distances
- Ordination: PCoA, NMDS, PCA
- Differential abundance: LEfSe, ANCOM-BC, DESeq2
- Visualization: stacked bar, heatmap, boxplot, phylogenetic tree
- R packages: `phyloseq`, `vegan`, `microbiome`, `mia`

### Talk 11 — Genome Analysis
- GWAS fundamentals: SNPs, LD, population structure
- Data formats: PLINK, VCF, BED
- GWAS workflow: quality control, association testing
- Visualization: Manhattan plot, QQ plot
- SNP annotation & functional prediction
- Linkage disequilibrium analysis
- Polygenic risk score (PRS) basics
- R packages: `SNPRelate`, `gaston`, `qqman`

### Talk 12 — Single-cell RNA-seq Analysis
- scRNA-seq technologies & data characteristics (sparsity, dropout)
- Seurat standard workflow:
  - Quality control & filtering
  - Normalization & scaling
  - Dimensionality reduction (PCA, UMAP, t-SNE)
  - Clustering & cluster marker identification
- Cell type annotation: automated & manual methods
- Differential expression between clusters/conditions
- Trajectory inference (Monocle / Slingshot)
- Data integration (multiple samples, batch correction)

### Talk 13 — Spatial Transcriptomics
- Spatial transcriptomics technologies (10x Visium, MERFISH, Slide-seq)
- Data structure: spatial coordinates + gene expression
- Loading & preprocessing spatial data (Seurat, Giotto)
- Spatial feature expression visualization
- Spatial domain/cluster identification
- Integration with scRNA-seq reference data
- Spatial niche analysis & cell-cell communication

### Talk 14 — Network Analysis
- Biological network types: co-expression, PPI, GRN
- WGCNA: weighted gene co-expression network analysis
  - Soft thresholding & scale-free topology
  - Module detection & module-trait association
  - Hub gene identification
- Protein-protein interaction networks (STRING database)
- Network topology analysis: degree, betweenness, clustering coefficient
- Cytoscape integration & network visualization
- (Brief) Gene regulatory network inference

### Talk 15 — Multi-omics Integration & Case Studies
- Multi-omics integration strategies & challenges
- Methods: MOFA, mixOmics (sPLS-DA, DIABLO)
- Cross-omics correlation analysis
- Case Study 1: Transcriptomics + Clinical data
  - DEGs → survival association → biomarker discovery
- Case Study 2: Microbiome + Metabolomics
  - Correlation network, pathway integration
- Reproducible research with Quarto
- Best practices: version control, renv, project organization

### Talk 16 — Other Omics & Course Wrap-up
- Proteomics: mass spectrometry data, protein quantification, differential analysis
- Metabolomics: LC-MS/GC-MS, metabolite identification, pathway mapping
- Epigenomics: DNA methylation (array, BS-seq), ATAC-seq, ChIP-seq
- Emerging omics: epitranscriptomics, lipidomics, glycomics
- Course summary: key skills & concepts revisited
- Further learning resources & career directions

---

## Additional Omics Topics (for future expansion)
- Drug repurposing & target prediction
- Pharmacogenomics
- Metagenomic binning & MAG analysis
- CRISPR screening data analysis
- Immuno-oncology (TCR/BCR repertoire)
- Long-read sequencing (Nanopore, PacBio) data analysis

---

*Last updated: 2026-05-27 — Refined Talk 01 outline with explicit 10-step structure and QMD-to-PDF verification exercises; added \footnotesize to all prompt texts*

