# GSE81608 human islet scRNA-seq analysis

Analysis of the GSE81608 dataset: single-cell RNA-seq of human pancreatic islet
cells from non-diabetic and type 2 diabetic (T2D) organ donors.

- **Study:** RNA Sequencing of Single Human Islet Cells Reveals Type 2 Diabetes Genes (PMID 27667665)
- **Technology:** Fluidigm C1 capture, SMARTer Ultra Low RNA kit, Nextera XT libraries, Illumina HiSeq 2500
- **Alignment:** GRCh37 (hg19), CLC Genomics Workbench 7.0
- **Values:** RPKM only (no raw counts are available for this dataset)
- **Size:** 1,600 cells, 39,851 genes, 18 donors (12 non-diabetic, 6 T2D)
- **Cell types (author annotation):** alpha, beta, delta, PP

## Directory layout

```
network/
├── README.md
├── methods.md                               # Proposed hdWGCNA methods (research proposal)
├── data/
│   ├── GSE81608_human_islets_rpkm.txt.gz   # GEO RPKM matrix (Entrez gene IDs x cells)
│   ├── GSE81608_series_matrix.txt.gz       # GEO series matrix (cell metadata)
│   ├── GSE81608_seurat.rds                 # Step 1: unfiltered Seurat object
│   ├── GSE81608_qc_filtered.rds            # Step 2: QC-filtered object
│   ├── GSE81608_umap.rds                   # Step 3: normalised, PCA, clusters, UMAP
│   ├── GSE81608_doublet_check.rds          # Step 4: Step 3 object + doublet flags
│   ├── GSE81608_harmony.rds                # Step 5: Harmony embedding, clusters, UMAPs
│   └── GSE81608_data_info.txt
├── scripts/
│   ├── read_data.R      + read_data.sh      # Step 1
│   ├── qc_filtering.R   + qc_filtering.sh   # Step 2
│   ├── umap.R           + umap.sh           # Step 3
│   ├── doublet_check.R  + doublet_check.sh  # Step 4
│   ├── harmony.R        + harmony.sh        # Step 5
│   └── logs/                                # SLURM logs (<step>_<jobid>.log)
└── results/
    ├── GSE81608_data_info.txt               # Step 1 text summary
    ├── GSE81608_qc_summary.pdf              # Step 1 QC overview plots
    ├── qc_filtering/                        # Step 2 outputs
    ├── umap/                                # Step 3 outputs
    ├── doublet_check/                       # Step 4 outputs
    └── harmony/                             # Step 5 outputs
```

## Running

Each step is an R script with a matching SLURM wrapper. Run them in order from
`scripts/`:

```bash
cd /home/alan.culligan/network/scripts
sbatch read_data.sh       # Step 1
sbatch qc_filtering.sh    # Step 2
sbatch umap.sh            # Step 3
sbatch doublet_check.sh   # Step 4
sbatch harmony.sh         # Step 5 (needs the Step 4 object)
```

All wrappers load `R/4.4.2-gfbf-2024a`, `R-bundle-CRAN/2024.11-foss-2024a` and
`R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2` (Seurat 5.1.0), and request
4 CPUs, 32 GB and 30 minutes (1 hour for `harmony.sh`). Each step takes
under 2 minutes. Paths and
parameters are set at the top of each R script.

## Donors

| Donor | Condition | Age | Sex | Ethnicity | Cells after QC |
|---|---|---|---|---|---|
| Non T2D 1 | non-diabetic | 23 | M | AA | 77 |
| Non T2D 2 | non-diabetic | 32 | F | C | 24 |
| Non T2D 3 | non-diabetic | 23 | F | C | 44 |
| Non T2D 4 | non-diabetic | 56 | F | C | 31 |
| Non T2D 5 | non-diabetic | 27 | M | AA | 57 |
| Non T2D 6 | non-diabetic | 68 | M | C | 50 |
| Non T2D 7 | non-diabetic | 29 | M | C | 29 |
| Non T2D 8 | non-diabetic | 60 | M | AI | 45 |
| Non T2D 9 | non-diabetic | 24 | F | C | 37 |
| Non T2D 10 | non-diabetic | 43 | M | C | 49 |
| Non T2D 11 | non-diabetic | 31 | F | H | 63 |
| Non T2D 12 | non-diabetic | 56 | M | AA | 144 |
| T2D 1 | T2D | 57 | M | AA | 90 |
| T2D 2 | T2D | 37 | F | C | 82 |
| T2D 3 | T2D | 55 | F | C | 69 |
| T2D 4 | T2D | 41 | F | H | 306 |
| T2D 5 | T2D | 42 | M | C | 175 |
| T2D 6 | T2D | 51 | M | H | 224 |

Ethnicity codes are as given in GEO (AA, AI, C, H). Cell numbers per donor are
very uneven (24–306): T2D 4, T2D 5, T2D 6 and Non T2D 12 together contribute
849 of the 1,596 cells.

## Step 1 – Read data (`read_data.R`)

Builds a Seurat object from the GEO files and writes a summary of the data.

1. Reads the RPKM matrix (genes x cells, first column Entrez gene ID).
2. Maps Entrez IDs to gene symbols with `org.Hs.eg.db`. 37,993 of 39,851 IDs
   map; the 1,858 unmapped IDs keep the Entrez ID as the feature name.
   Underscores are replaced with dashes (Seurat requirement). The mapping is
   stored in the RNA assay feature metadata (`entrez_id`, `symbol`, `feature`).
3. Reads cell metadata from the series matrix: `donor_id`, `condition`, `age`,
   `ethnicity`, `gender`, `cell_type`, `tissue`, `geo_accession`.
4. Creates the Seurat object with RPKM values in the `counts` layer (no cells or
   genes removed) and adds QC metrics: `percent.mt` (37 chrMT genes),
   `log10_nCount_RNA`, `log10_nFeature_RNA`, `log10GenesPerCount` and
   `pct_counts_top20`.

Because the values are RPKM, `nCount_RNA` is summed RPKM per cell, not
library size.

**Outputs:** `data/GSE81608_seurat.rds`, `results/GSE81608_data_info.txt`,
`results/GSE81608_qc_summary.pdf`.

**Data summary (before QC):**
- Genes detected per cell: median 5,688 (2,406–11,142).
- Summed RPKM per cell: median 538,526.
- Mitochondrial RPKM: median 16.6% (1.9–45.4%).
- Top genes by mean RPKM are GCG, TTR, INS, mitochondrial transcripts, PPY and SST, as expected for islet cells.

## Step 2 – QC filtering (`qc_filtering.R`)

Per-donor, MAD-based filtering, following the APAP QC pipeline
(`pipelines/qc_filtering`: `bin/qc_module.R` and `modules/qc_per_sample.nf`).
Each donor is treated as a sample.

**Filters:** two-sided median ± 5 MAD (MAD constant 1.4826), computed separately
for each donor, on:
- `nFeature_RNA`
- `log1p_total_counts` (log1p of summed RPKM)
- `log1p_n_genes_by_counts`
- `pct_counts_in_top_20_genes`

Not used for filtering:
- **Mitochondrial content:** it is only plotted, because it differs strongly between islet cell types.
- **Hard `nFeature_RNA` bounds:** these are switched off.

**Result:** 1,596 of 1,600 cells kept (99.8%).
- 4 cells were removed: 1 failed `log1p_total_counts` and 3 failed `log1p_n_genes_by_counts`.
- The removed cells came from Non T2D 9, T2D 1, T2D 2 and T2D 5, one each.
- By cell type, the removed cells were 2 alpha, 1 beta and 1 PP.

The data were already clean. GEO holds 1,600 cells; the published analysis
reports 1,492 alpha, beta, delta and PP cells.

**Outputs (`results/qc_filtering/`):**
- `per_donor/<donor>/`: histograms with the donor's MAD bounds, a thresholds report and a log.
- `per_condition/<condition>/`: pooled histograms, post-filter violin, ridge and scatter plots, and a thresholds report.
- `QC_plots_combined.pdf`, `Thresholding_combined.txt` and `seurat_processing_log.txt`.

The filtered object is `data/GSE81608_qc_filtered.rds`. Its filter settings
and per-donor bounds are stored in `Misc(seu, "qc_filtering")`.

## Step 3 – PCA, PC selection and UMAP (`umap.R`)

**Method:**
1. **Normalisation:** RPKM is already corrected for gene length and sequencing
   depth, so the `data` layer is set to `log1p(RPKM)`. `LogNormalize` is not
   used, because it would rescale every cell to the same total.
2. **Variable features:** the top 2,000 genes (vst).
   - Islet hormones: SST, PPY, NPY and IAPP.
   - Exocrine (acinar and ductal) genes: REG1A, SPINK1, LCN2 and MMP7.
   - Stromal genes: COL1A1, SPARC and LUM.
   - Inflammatory and stress genes: CCL2, CXCL1, CXCL8 and MT2A.
   - The non-islet genes may reflect ambient RNA or contamination from cells
     around the islets.
3. **Scaling and PCA:** the data are scaled and 50 PCs are computed (seed 42).
4. **PC selection from the elbow:**
   - (a) Cumulative rule: the first PC at which the PCs together hold more than
     90% of the standard deviation and the PC itself holds less than 5%.
   - (b) Drop rule: the last PC where the drop to the next PC is still more
     than 0.1%.
   - The estimate is min(a, b). These two rules are the HBC training heuristics.
   - (c) A geometric knee of the elbow curve is also reported for comparison.
   - To override the estimate, set `n_pcs` at the top of the script.
5. **Clustering and UMAP:** `FindNeighbors` and `FindClusters` (Louvain,
   resolution 0.5), then `RunUMAP` (uwot, 30 neighbours, seed 42), all on the
   chosen PCs.
6. **Sensitivity check:** UMAPs are also computed at 5, 10, 15, 20 and 30 PCs.

**PC selection result:**

| Estimate | PCs |
|---|---|
| (a) Cumulative | 44 |
| (b) Drop | 10 |
| (c) Knee | 7 |
| **Used, min(a, b)** | **10** |

The first 10 PCs explain 10.5% of the variance of the scaled data. PC1–PC4
carry most of the structure, and the curve is flat after about PC10. The
sensitivity UMAPs show the same main structure from 10 to 30 PCs. With 5 PCs,
the delta and PP cells are not fully separated.

**UMAP result:**
- **Cell type drives the structure.** Alpha, beta, delta and PP cells form
  clearly separate groups, and the islet markers confirm this: GCG in alpha,
  INS in beta, SST in delta and PPY in PP cells.
- **Gender, age and ethnicity do not form their own clusters.**
- **Condition is mixed within each cell-type group.** Within the alpha and beta
  groups there are moderate differences in the proportion of T2D cells between
  subclusters, which is expected.
- **10 clusters at resolution 0.5:**

| Cluster | Cells | Main cell type | T2D share |
|---|---|---|---|
| 0 | 403 | alpha | 76% |
| 1 | 312 | alpha | 51% |
| 2 | 290 | beta | 78% |
| 3 | 175 | beta | 29% |
| 4 | 160 | alpha | 44% |
| 5 | 81 | PP | 62% |
| 6 | 55 | mixed (30 alpha, 15 beta) | 58% |
| 7 | 51 | mixed (30 alpha, 14 beta) | 61% |
| 8 | 47 | delta | 47% |
| 9 | 22 | alpha (18 of 22 from Non T2D 1) | 5% |

**Points to follow up:**
- **The alpha and beta subclusters split by donor.**
  - Clusters 0 (alpha) and 2 (beta) contain 593 of their 693 cells from four
    donors: T2D 4, T2D 5, T2D 6 and Non T2D 12. Clusters 1, 3 and 4 hold the
    other donors.
  - The T2D enrichment in clusters 0 and 2 therefore mostly reflects these four
    donors. Non T2D 12 groups with the T2D donors, and T2D 1–3 group with the
    non-diabetic donors.
  - Cells in clusters 0 and 2 also detect fewer genes (median 5,051 vs 6,156).
  - This may be a donor or processing-batch effect rather than a disease
    effect. Step 5 tests this with Harmony.
- **Clusters 6 and 7 contain a mix of alpha and beta cells.** Step 4 shows
  they are endocrine cells combined with ductal (cluster 6) or stellate and
  endothelial (cluster 7) cells, not alpha + beta doublets.
- **Cluster 9 is almost entirely one donor (Non T2D 1).** It merges into the
  main alpha cluster after Harmony (Step 5).

**Outputs (`results/umap/`):**
- `elbow_plot.pdf`:
  - Variable feature plot.
  - Elbow plot with the PC estimates marked.
  - Variance explained per PC and cumulatively.
  - PC heatmaps for PC1–PC24.
  - PC1 vs PC2 and PC3 vs PC4, coloured by donor, condition and cell type.
- `umap_plots.pdf`:
  - Overview by donor, condition, cell type and cluster, then full-page UMAPs by donor, condition and cell type.
  - Cell type split by condition and by donor.
  - Gender, ethnicity and age.
  - QC metrics.
  - Marker genes.
  - Sensitivity to the number of PCs.
- `pc_selection.txt`: the PC estimates and the variance table for each PC.
- `cluster_tables.txt`: clusters by cell type, by condition and by donor, and donors by cell type.

The processed object is `data/GSE81608_umap.rds`. It has the `pca` and `umap`
reductions, `seurat_clusters`, and the settings in `Misc(seu, "umap")`.

## Step 4 – Doublet check (`doublet_check.R`)

Fluidigm C1 capture sites can hold two cells. DoubletFinder and scDblFinder
expect raw UMI counts, and this dataset has only RPKM, so the script uses three
RPKM-compatible checks.

**1. Hormone co-expression.**
- For each cell, the highest off-type hormone (GCG, INS, SST, PPY) is divided
  by the hormone of its annotated type.
- Ambient hormone is around 0.1–1% of the own hormone. A heterotypic doublet
  is around 10–50%.
- Flag: ratio > 0.1, or a dominant hormone that does not match `cell_type`.

**2. Simulated doublets.**
- 399 artificial heterotypic doublets (pN = 0.25) are made by averaging the
  RPKM profiles of random cell pairs of different annotated types.
- They are projected into the Step 3 PCA space. The projection reproduces the
  stored PCs exactly (correlation 1.000).
- Each cell is scored by the fraction of artificial doublets among its 20
  nearest neighbours.
- Flag: score > 0.4, i.e. 2 x the overall artificial share of 0.2.

**3. Non-endocrine signal.**
- Each cell gets the mean log1p(RPKM) of five marker sets:
  - acinar: PRSS1, CPA1, CELA3A, CTRB1, REG1A
  - ductal: KRT19, CFTR, SPP1
  - stellate: COL1A1, SPARC, PDGFRB
  - endothelial: PECAM1, PLVAP
  - immune: PTPRC
- Flag: any set > 3. At this level the sets separate cleanly, and the
  widespread low-level acinar signal is not flagged.

**Combined call (`doublet_call`), in order of priority:**
1. endocrine + non-endocrine (check 3)
2. endocrine + endocrine (checks 1 and 2 agree)
3. hormone only
4. simulation only
5. singlet

**Result:**

| Call | Cells |
|---|---|
| singlet | 1,428 |
| endocrine + non-endocrine | 107 |
| simulation only | 45 |
| hormone only | 10 |
| endocrine + endocrine | 6 |

- **Endocrine + endocrine doublets are rare.**
  - Only 18 cells (1.1%) co-express a second hormone above 10% of their own,
    almost all beta + GCG (13). No cell has a dominant hormone that disagrees
    with its label.
  - Only 6 cells are supported by both hormone and simulation evidence.
  - The 45 simulation-only cells sit mostly in the small delta and PP groups,
    which lie close to where artificial alpha + delta and alpha + PP doublets
    land. Treat these as low confidence.
- **Endocrine + non-endocrine is the main problem.**
  - 107 cells are flagged: 48 ductal, 36 stellate, 11 endothelial, 10 acinar
    and 2 immune.
  - These include 85% of cluster 6 and 90% of cluster 7.
  - These cells keep their own hormone at about half the level of pure cells
    (median GCG about 107,000–141,000 RPKM vs 189,000–235,000), express a full
    ductal or stellate programme, and detect more genes. This is the expected
    profile of a doublet of an endocrine cell and a non-endocrine cell.
- **Acinar transcripts are widespread at a low level.** CELA3A, CPA1, PRSS1
  and REG1A appear across alpha cluster 4 and beta cluster 3 (see Step 5
  markers), but mostly below the flag threshold. This looks like ambient
  exocrine RNA in the islet preparation rather than doublets.

**Outputs (`results/doublet_check/`):**
- `doublet_check.pdf`:
  - Hormone scatter plots and the ratio distribution.
  - Artificial doublets in PCA space and the score distribution.
  - UMAPs of the calls.
  - Proportion of calls by cluster, cell type, donor and condition.
  - Marker dot plots by cluster and by cell type.
  - Non-endocrine score violins.
- `doublet_summary.txt`: the thresholds, the counts, and the proportions by cluster, cell type, donor and condition.
- `doublet_flags.csv`: one row per cell with its values and calls.

## Step 5 – Harmony integration and comparison (`harmony.R`)

**Method:**
1. **Harmony:** run on the 10 Step 3 PCs, grouping by `donor_id` (theta 2).
2. **Clustering and UMAP:** done for the unintegrated PCA and for the Harmony
   embedding with identical settings (10 PCs, resolution 0.5, 30 neighbours,
   seed 42), plus a resolution grid of 0.1–1 for the clustering trees.
3. **Doublets kept:** all cells are kept by default. To remove doublet calls
   first, set `exclude_calls`, for example to
   `c("endocrine + endocrine", "endocrine + non-endocrine")`.

Harmony changes only the embedding, not the expression values. Differential
expression should still use uncorrected data with donor taken into account
(for example pseudobulk per donor). Each donor belongs to one condition, so
integrating by donor can also remove part of a real condition difference from
the embedding.

**Metrics:**

| Metric | Unintegrated | Harmony | Better if |
|---|---|---|---|
| Clusters (res 0.5) | 10 | 8 | – |
| Median donor diversity (kNN, k = 30; 10.59 if fully mixed) | 4.25 | 6.62 | higher |
| Median condition diversity (1.93 if fully mixed) | 1.72 | 1.92 | higher |
| Median cell type diversity | 1.00 | 1.00 | stays at 1 |
| ASW cell_type | 0.325 | 0.365 | higher |
| ASW donor_id | -0.147 | -0.196 | lower |
| ARI clusters vs cell_type | 0.362 | 0.491 | higher |
| NMI clusters vs cell_type | 0.408 | 0.487 | higher |

Harmony improved donor mixing without merging cell types. The four cell types
stay fully separate, and the clusters match the annotation better.

**Clusters before and after:**

| Unintegrated | Harmony | Interpretation |
|---|---|---|
| 2 + 3 (beta) | 1 (463 beta) | Beta split was donor-driven; merges into one cluster |
| 1 + most of 9 (alpha) | 0 (491 alpha) | Main alpha cluster |
| 0 (alpha) | 2 (258) and part of 0 | Partly persists, weak markers (see below) |
| 4 (alpha) | 3 (124 alpha) | Alpha with acinar transcripts |
| 6 | 6 (60) | Endocrine + ductal |
| 7 | 5 (63) | Endocrine + stellate / endothelial |
| 5 | 4 (84) | PP |
| 8 | 7 (47) | Delta |

**Why there are several clusters within one annotated cell type:**

The `cell_type` labels are the authors' identities. Seurat clusters are
unsupervised partitions at an arbitrary resolution (0.5). They split a cell
type wherever there is enough structure, whether biological, technical or
donor-related. The marker comparison within alpha and beta
(`subcluster_markers.csv`) shows what drives each split:

| Subclusters | Top markers | Driver |
|---|---|---|
| Unintegrated beta 2 vs 3 | Beta 2: RPL21 pseudogenes, LINC00486. Beta 3: PRSS1, CPA2, PNLIP, REG1A | Donor / technical; acinar contamination in beta 3. Merged by Harmony |
| Unintegrated alpha 0 vs 1 | Alpha 0: RPL21 pseudogenes, LINC00486 (log2FC ≤ 1). Alpha 1: weak | Donor / technical |
| Harmony alpha 0 vs 2 | Weak (log2FC ≤ 1.8, no known alpha-cell genes) | Residual donor / technical |
| Unintegrated alpha 4 = Harmony alpha 3 | CELA3A/B, CPA1, CTRB1/2, CLPS, PNLIP | Acinar contamination |
| Unintegrated 6 = Harmony 6 | MMP7, CLDN2, CEACAM6/7, FUT3 | Ductal doublets |
| Unintegrated 7 = Harmony 5 | POSTN, COL6A3, SPARC, VCAN, FOXF2 | Stellate / endothelial doublets |

About the donor-driven subclusters:
- The donor-driven groups are marked by pseudogenes (RPL21P*) and LINC00486,
  together with fewer genes detected. This points to a library or processing
  batch shared by T2D 4–6 and Non T2D 12, rather than biology.
- Harmony alpha 2 still holds 73% of its cells from these four donors.

So far, none of the alpha or beta subclusters shows a distinct biological
alpha- or beta-cell state. They are explained by donor or batch,
contamination, or doublets.

**Outputs (`results/harmony/`):**
- `harmony_comparison.pdf`, comparing unintegrated and Harmony:
  - UMAPs by donor, condition, cell type, clusters and doublet call.
  - Cell type split by condition and by donor.
  - Donor and cell-type composition of each cluster.
  - Integration metrics.
  - Clustering trees.
  - Dot plots of the subcluster markers.
- `harmony_summary.txt`: the settings, metrics and cluster tables, old vs new clusters, and the top 10 markers per subcluster.
- `subcluster_markers.csv`: all significant subcluster markers (adjusted p < 0.05), before and after.

The processed object is `data/GSE81608_harmony.rds`. It has the reductions
`pca`, `harmony`, `umap_pca` and `umap_harmony`, the columns `clusters_pca`
and `clusters_harmony` (`seurat_clusters` = Harmony clusters), the
resolution-grid columns, and the settings in `Misc(seu, "harmony")`.

## Planned – hdWGCNA co-expression network analysis of beta cells

**Aim:** use hdWGCNA (Morabito et al. 2023) to find gene co-expression modules
in beta cells and compare them between T2D and non-diabetic donors.

**Input:** the 502 beta cells in the QC-filtered object (207 from 12
non-diabetic donors, 295 from 6 T2D donors), on log1p(RPKM).

**Planned comparisons:**
- **Joint beta-cell network** built from all donors, then differential module
  eigengene testing between T2D and non-diabetic donors at the donor level, so
  that cells from the same donor are not counted as independent replicates.
- **Module preservation:** build the network in non-diabetic beta cells and
  test which modules are preserved, weakened or lost in T2D beta cells, and the
  reverse.
- **Consensus network** across the two conditions, plus differential
  connectivity (kME) of hub genes between conditions.
- **Module–trait relationships** with condition, age, sex and technical
  covariates (genes detected, mitochondrial %, acinar signal).
- **Cell-type specificity:** project the beta-cell modules into alpha, delta
  and PP cells.
- **Functional annotation:** GO/KEGG/Reactome enrichment of each module, and
  overlap with beta-cell identity, stress and T2D GWAS gene sets.

**Points to handle:**
- Values are RPKM, so metacells are averaged RPKM followed by log1p, not
  hdWGCNA's default `NormalizeMetacells` (LogNormalize).
- Donor is nested in condition. Harmony on `donor_id` can remove real condition
  signal, so harmonised eigengenes (hMEs) are used only for visualisation, and
  condition tests use non-harmonised MEs with donor accounted for.
- The likely processing batch (T2D 4–6 and Non T2D 12, see Steps 3 and 5)
  overlaps with condition. Results are repeated without these donors.
- Cell numbers per donor are uneven (beta: 7–144), so each donor's metacell
  contribution is capped.

hdWGCNA (and GeneOverlap) are not yet installed in the R 4.4.2 module
environment. WGCNA 1.73, UCell 2.10.1, harmony 1.2.1, lme4 1.1.35.5 and
clusterProfiler 4.14.3 are available.

The full proposed methods are in `methods.md`.

## Notes

- All values are RPKM. Methods that need raw counts (for example SCTransform,
  DoubletFinder's count-based steps, or count-based differential expression
  such as DESeq2 or edgeR on single cells) are not directly applicable. Use
  log-RPKM-based or pseudobulk approaches.
- Mitochondrial content is high in this dataset (median 16.6%) and differs
  between cell types, so it is not used as a filter.
- Cell type labels come from the original authors (GEO `cell subtype`).
