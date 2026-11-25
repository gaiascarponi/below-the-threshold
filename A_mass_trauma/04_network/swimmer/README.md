# SWIMmeR: switch genes, Axis A (WTC)

## Thesis run: `seed_18/`

SWIMmeR (Paci et al., 2017) was run in R 4.6.1 with `set.seed(18)`, set in `main.R` and before the `kmeans()` call in `getClustering.R`, so that the k-means clustering is reproducible. The SWIMmeR source code is not included in this repository.

The result was also compared with an earlier run on the same input (expression matrix of 24,547 genes x 324 samples, group lists of 201 and 123 samples) without a fixed seed: the 51 switch genes, the network and the cluster sizes are identical, so the result does not depend on the k-means initialization.

## Parameters

- SWIMmeR internal filter: unpaired t-test, FDR correction, adjusted p < 0.05, fold change > 1.1, IQR filter 20%, at most 75% zeros
- Correlation: Pearson, threshold at the 0.8 percentile (r >= 0.3452), adjusted p of the correlation < 0.05
- k-means: k = 3, iter.max = 100, 5 repetitions

## Result

Network of 182 nodes and 6,878 edges, clusters of 72, 72 and 38 nodes, k-means partition score 68.6. **51 switch genes** (region R4, APCC < 0), 50 of which are up-regulated in PTSD.

## Main files (`seed_18/switch/txtFile/`)

- `switch.txt`: the 51 switch genes
- `attribute_switch.txt`: attributes of the switch genes (degree, APCC, z, P, cluster, logFC, adjusted p)
- `attribute.txt`: attributes of all 182 nodes of the network
- `CorrelationNetwork.txt`, `CartographyNetwork.txt`: network edges
- `nn_pos_switch.txt`, `nn_neg_switch.txt`: positively and negatively correlated neighbors of the switch genes
- Figures in `seed_18/switch/figure/`: heat cartography map, APCC distribution, switch gene heatmap, connected components
