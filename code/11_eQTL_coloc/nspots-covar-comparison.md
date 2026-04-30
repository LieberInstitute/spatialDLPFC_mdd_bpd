# nspots Covariate Comparison

This note compares the baseline tensorQTL summaries against the nspots-adjusted tensorQTL summaries generated from `processed-data/11_eQTL_coloc/nspots`.

Bottom line: adding `nspots` had a modest but systematic effect, mostly reducing discovered eGenes/signals, especially in sex-stratified groups. It did not change any DEG+GWAS overlap counts in these summaries.

## Independent eQTL Summary

These are the side-by-side columns for filtered independent signals.

| group | context | eGenes_base | eGenes_nspots | delta | DEG_base | DEG_nspots | delta_DEG | DEG_GWAS_base | DEG_GWAS_nspots |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all | Astro | 586 | 578 | -8 | 32 | 35 | +3 | 1 | 1 |
| all | Inhb | 271 | 272 | +1 | 23 | 23 | 0 | 0 | 0 |
| all | L2.3 | 994 | 960 | -34 | 35 | 34 | -1 | 1 | 1 |
| all | L4 | 215 | 216 | +1 | 20 | 20 | 0 | 1 | 1 |
| all | L5 | 807 | 799 | -8 | 30 | 32 | +2 | 1 | 1 |
| all | L6 | 491 | 477 | -14 | 31 | 30 | -1 | 0 | 0 |
| all | Micro.Vasc | 23 | 24 | +1 | 3 | 3 | 0 | 0 | 0 |
| all | Oligo | 85 | 79 | -6 | 11 | 12 | +1 | 0 | 0 |
| male | Astro | 104 | 93 | -11 | 9 | 7 | -2 | 0 | 0 |
| male | Inhb | 43 | 43 | 0 | 5 | 4 | -1 | 0 | 0 |
| male | L2.3 | 140 | 129 | -11 | 8 | 7 | -1 | 1 | 1 |
| male | L4 | 24 | 29 | +5 | 2 | 3 | +1 | 0 | 0 |
| male | L5 | 126 | 115 | -11 | 7 | 6 | -1 | 0 | 0 |
| male | L6 | 92 | 81 | -11 | 5 | 4 | -1 | 0 | 0 |
| male | Micro.Vasc | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| male | Oligo | 4 | 4 | 0 | 1 | 1 | 0 | 0 | 0 |
| female | Astro | 109 | 104 | -5 | 10 | 10 | 0 | 0 | 0 |
| female | Inhb | 49 | 47 | -2 | 7 | 5 | -2 | 0 | 0 |
| female | L2.3 | 176 | 164 | -12 | 12 | 13 | +1 | 0 | 0 |
| female | L4 | 45 | 39 | -6 | 3 | 3 | 0 | 0 | 0 |
| female | L5 | 149 | 134 | -15 | 12 | 10 | -2 | 0 | 0 |
| female | L6 | 67 | 61 | -6 | 7 | 7 | 0 | 0 | 0 |
| female | Micro.Vasc | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| female | Oligo | 2 | 2 | 0 | 1 | 1 | 0 | 0 | 0 |

## Group-Level Effect

| result | group | base | nspots | delta | relative |
|---|---:|---:|---:|---:|---:|
| map_cis eGenes | all | 3359 | 3296 | -63 | -1.9% |
| map_cis eGenes | male | 533 | 495 | -38 | -7.1% |
| map_cis eGenes | female | 593 | 547 | -46 | -7.8% |
| independent signals | all | 3472 | 3405 | -67 | -1.9% |
| independent signals | male | 535 | 496 | -39 | -7.3% |
| independent signals | female | 598 | 552 | -46 | -7.7% |

## DEG-View Overlaps

Aggregated across contexts within each group, overlap changes were small. These counts are not mutually exclusive across DEG views.

| group | DEG view | overlap_base | overlap_nspots | delta |
|---|---:|---:|---:|---:|
| all | broad_interaction | 241 | 242 | +1 |
| all | context_localized | 26 | 28 | +2 |
| all | context_and_sex_specific:female | 23 | 25 | +2 |
| all | context_and_sex_specific:male | 6 | 6 | 0 |
| all | sex_specific:female | 89 | 91 | +2 |
| all | sex_specific:male | 62 | 63 | +1 |
| male | broad_interaction | 46 | 41 | -5 |
| male | context_localized | 8 | 7 | -1 |
| male | sex_specific:male | 10 | 8 | -2 |
| female | broad_interaction | 62 | 58 | -4 |
| female | context_localized | 6 | 7 | +1 |
| female | context_and_sex_specific:female | 5 | 6 | +1 |
| female | sex_specific:female | 8 | 8 | 0 |

The biggest individual DEG-view overlap shifts were small: Astro/all `sex_specific:female` increased from 16 to 20, Astro/all `broad_interaction` increased from 40 to 43, and several sex-stratified rows moved by only 1-2 genes.

## Interpretation

Adding `nspots` is not a dramatic re-wiring of the eQTL results. It mostly trims signal counts, with the largest relative effect in male/female-only runs, around 7-8% fewer eGenes/signals overall. The all-donor runs changed much less, around 2% fewer.

The DEG overlap story is even more stable. DEG-view overlaps move by only a few genes per view, and `n_DEG_GWAS` is unchanged everywhere. So `nspots` matters statistically enough to reduce some marginal discoveries, especially in smaller sex-stratified models, but it does not materially change the biological overlap conclusions from these summaries.
