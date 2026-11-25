# Dictionary: references for gene identifier harmonization and gene sets

Used by `Transferability/00_harmonize_gene_ids.R` and `Transferability/01_signature_transfer_gsea.R` (functions `build_gene_dictionary()`, `to_entrez()` and `harmonize_deg_table()` in `00_functions.R`) to match genes across axes through **NCBI Entrez Gene IDs** instead of gene symbols, and by the Level 3 scripts for gene sets and cell-type markers.

| File | Source | Version |
|---|---|---|
| `hgnc_complete_set.txt` | HGNC complete set, https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | downloaded on 8 October 2026 |
| `Homo_sapiens.gene_info.gz` | NCBI Gene FTP, https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | NCBI update of 8 October 2026, downloaded the same day |
| `msigdb_H.rds`, `msigdb_C5_GO_BP.rds` | MSigDB Hallmark and GO Biological Process gene sets, used in `03_leading_edge.R` | obtained with the `msigdbr` R package (version in `sessionInfo_R.txt`) |
| `MCPcounter_genes.txt` | marker genes of the MCP-counter package (v1.2.0, GPL-3), used in `04_cell_composition.R` | blood cell populations only |

## Conversion rules (in order)

1. Ensembl ID (Axis C) to Entrez, unique matches only;
2. official NCBI symbol to Entrez;
3. aliases and previous symbols (NCBI synonyms, HGNC `alias_symbol` and `prev_symbol`) to Entrez, only if the alias maps to a single gene (otherwise discarded as ambiguous);
4. `LOC<number>` to the number, if it is an existing Entrez ID;
5. gene names corrupted by Excel (e.g. `1-Mar`, `2-Sep`) are excluded.

Replacing these files with newer versions may change the results.
