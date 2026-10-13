### Script to download the GSE221156 'reintegration of betacells dataset' from
### CellxGene - which includes authors processed data (QC, normalised counts, clusterings and
### additional metadata)

# Path configuration
Root <- '/Users/alan/Desktop/Bioinformatics/net-hdWGCNA_betacells'
Out_dir <- file.path(Root, '1.Data', 'GSE221156')
Out_file <- file.path(Out_dir, 'GSE221156_beta_cellxgene.h5ad')
CellxGene_link <- 'https://datasets.cellxgene.cziscience.com/faacd98e-af54-459d-8a7a-74427766bc81.h5ad'

# Curl dataset
library(curl)
dir.create(Out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(Out_file)) {
  curl_download(CellxGene_link, Out_file, quiet = FALSE)
}
