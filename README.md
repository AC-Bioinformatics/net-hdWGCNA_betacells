# hdWGCNA co-expression network analysis of human beta cells

Gene co-expression network analysis (hdWGCNA) of human pancreatic beta cells
from the [GSE221156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156) single-cell RNA-seq dataset, comparing donors who are
non-diabetic (ND), prediabetic (PD) and type 2 diabetic (T2D).

## Dataset

- **GEO accession:** GSE221156 (Bandesh et al.)
- **Source file:** the beta-cell subset downloaded from CELLxGENE as an `.h5ad`
- **Technology:** 10x Genomics 3' (v2 and v3 chemistry)
- **Counts:** raw UMI counts from `raw.X`. The authors' normalised `X` is not used.
- **Cells:** beta cells only, with the authors' beta-cell subclusters (`Clusters`)
  and the UMAP from their beta-cell reintegration

## Directory layout

```
network/
├── README.md
├── 1.Data/
│   └── GSE221156/
│       ├── GSE221156_beta_cellxgene.h5ad   # CELLxGENE download (not on GitHub)
│       └── GSE221156_beta_seurat.rds       # Output of step 1 (not on GitHub)
├── 2.Scripts/
│   ├── 1.read_GSE221156.R                  # Step 1: h5ad -> Seurat
│   └── 2.loadRDS.R                         # Step 2: inspect Seurat metadata
├── 3.Results/
└── data/                                   # Earlier GSE81608 work, not used here
```

## Workflow

### Step 1: read the h5ad into Seurat (`1.read_GSE221156.R`)

- Reads the raw UMI counts from `raw/X` and checks that they are integers
- Uses gene symbols from `var$feature_name` (Ensembl ID if no symbol), replaces
  underscores with dashes and makes duplicated symbols unique
- Keeps all cell metadata except ontology IDs, constant columns and the authors'
  `nCount_RNA`/`nFeature_RNA`, which Seurat recomputes from the same counts
- Adds short labels used downstream:

  | Column      | Values / source                                    |
  |-------------|----------------------------------------------------|
  | `condition` | `ND`, `PD`, `T2D` (from `disease`)                 |
  | `donor`     | islet donor (from `Islet`)                         |
  | `chemistry` | `v2`, `v3` (from `assay`)                          |
  | `age`       | years, integer (from `development_stage`)          |

- Log-normalises the counts (`NormalizeData`)
- Stores the authors' UMAP as the `umap_authors` reduction
- Prints cells per condition, donors per condition, beta cells per donor,
  chemistry by condition and author subclusters by condition
- Saves `1.Data/GSE221156/GSE221156_beta_seurat.rds`

### Step 2: inspect metadata (`2.loadRDS.R`)

Loads the Seurat object and writes the metadata columns and structure to
`GSE221156_summary.txt`.

### Next steps

hdWGCNA network construction on the beta cells, comparing modules across
ND, PD and T2D.

## Data availability

Data objects (`.h5ad`, `.rds`) and downloaded raw files are kept on the TGX
HPC and are excluded by `.gitignore`, because GitHub rejects files over 100 MB.
The h5ad can be downloaded from CELLxGENE (GSE221156 beta-cell subset).

## Requirements

R with `Seurat`, `rhdf5` and `Matrix`. `hdWGCNA` is needed for the network
analysis.
