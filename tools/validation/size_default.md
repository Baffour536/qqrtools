# Size of pointwise tests with the default bandwidth and degree

Null designs (no effect; N: normal X, S: skewed X), 300 replications, 5 x 5 cells, h = "slope", degree = 2, nominal level 5%; bootstrap with R = 100. The last two columns give the share of fitted surfaces with at least one significant cell before and after `qqr_padjust(method = "BH")`.

| design | n | se_type | pooled_rejection_rate | max_cell_rate | any_cell_significant | any_cell_significant_BH |
|---|---|---|---|---|---|---|
| N | 250 | nid | 0.093 | 0.170 | 0.820 | 0.407 |
| N | 250 | ker | 0.035 | 0.090 | 0.483 | 0.060 |
| N | 250 | boot | 0.050 | 0.077 | 0.577 | 0.087 |
| N | 1000 | nid | 0.079 | 0.153 | 0.747 | 0.280 |
| N | 1000 | ker | 0.042 | 0.080 | 0.513 | 0.080 |
| N | 1000 | boot | 0.050 | 0.080 | 0.567 | 0.087 |
| S | 250 | nid | 0.076 | 0.157 | 0.740 | 0.263 |
| S | 250 | ker | 0.042 | 0.093 | 0.517 | 0.083 |
| S | 250 | boot | 0.050 | 0.067 | 0.500 | 0.097 |
| S | 1000 | nid | 0.074 | 0.150 | 0.767 | 0.237 |
| S | 1000 | ker | 0.044 | 0.100 | 0.563 | 0.097 |
| S | 1000 | boot | 0.052 | 0.073 | 0.613 | 0.103 |
