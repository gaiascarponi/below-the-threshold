# Axis A: network and regulatory analysis (Level 2)

Three analyses of the Axis A signature (525 DEGs):

1. **Switch genes** with SWIMmeR. Input: normalized expression matrix of the 324 samples and the PTSD Yes/No sample lists (`SWIMmeR/project/PTSD_WTC/dataset/PTSD/`; the matrix is not included in this repository because of its size). Results in `swimmer/`.
2. **miRNA regulation** with MIENTURNET (web tool), on the up-regulated DEGs, the down-regulated DEGs and the switch genes. Exported tables in `mienturnet/`.
3. **Protein-protein interaction network** with STRING (web tool), on the DEGs and the switch genes. Exported networks and node tables in `string/`.

`input/` contains the gene lists submitted to the web tools: one gene symbol per line, ordered by adjusted p-value, written by `export_deg_lists()` in `A_mass_trauma/01_DEG.R` (`A_DEG_all.txt` 525 genes, `A_DEG_up.txt` 211, `A_DEG_down.txt` 314).
