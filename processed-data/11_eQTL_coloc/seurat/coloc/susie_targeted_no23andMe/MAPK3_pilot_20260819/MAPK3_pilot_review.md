# Reviewed MAPK3 SuSiE pilot

- Completed comparisons: 24 (8 contexts x 3 disorders).
- Primary fits converged: 48/48.
- Direct susie_rss and coloc::runsusie fits were numerically identical for every trait.
- Maximum estimate_s_rss: eQTL 0.000457812; GWAS 0.0423377.
- Maximum GWAS finite-reference effective-rank/B: 0.0396; no reliability flags.
- Naive n-row individual SuSiE was not used for scaling because residualization did not preserve covariate degrees of freedom.
- Residual-space compressed individual SuSiE and rank-aware RSS differed by at most 0.0413 PIP in the representative MDD fits.
- Stored tensorQTL RSS and rank-aware reconstructed RSS differed by at most 0.0944 PIP; all credible-set counts agreed.
- Primary later-target method: runsusie on stored dense tensorQTL statistics with rank-aware covariate-adjusted donor LD; trait-specific prefit coloc.susie.
- GWAS uses approved no-23andMe/current summary statistics and 503-EUR finite-reference correction.
- Raw donor LD remains a sensitivity; no nearPD, shrinkage, significance pruning, or reduced coverage was used.

## Primary MAPK3 results

    disorder context               hit1               hit2        PP3
      <char>  <char>             <char>             <char>      <num>
 1:       BD   astro chr16:30123335:T:C chr16:29946895:G:A 0.07871683
 2:      MDD   astro chr16:30123335:T:C chr16:30406798:G:A 0.91947870
 3:     SCZD   astro chr16:30123335:T:C chr16:29983601:C:T 0.99237203
 4:       BD    inhb chr16:30311847:G:C chr16:29946895:G:A 0.98546728
 5:      MDD    inhb chr16:30311847:G:C chr16:30406798:G:A 0.01645991
 6:     SCZD    inhb chr16:30311847:G:C chr16:29983601:C:T 0.99984464
 7:       BD    l2-3 chr16:30123335:T:C chr16:29946895:G:A 0.17545145
 8:      MDD    l2-3 chr16:30123335:T:C chr16:30406798:G:A 0.94097419
 9:     SCZD    l2-3 chr16:30123335:T:C chr16:29983601:C:T 0.88626922
10:       BD      l4 chr16:30123335:T:C chr16:29946895:G:A 0.16698477
11:      MDD      l4 chr16:30123335:T:C chr16:30406798:G:A 0.91515693
12:     SCZD      l4 chr16:30123335:T:C chr16:29983601:C:T 0.87444108
13:       BD      l5 chr16:30123335:T:C chr16:29946895:G:A 0.12674260
14:      MDD      l5 chr16:30123335:T:C chr16:30406798:G:A 0.88225114
15:     SCZD      l5 chr16:30123335:T:C chr16:29983601:C:T 0.88517600
16:       BD      l6                                               NA
17:      MDD      l6                                               NA
18:     SCZD      l6                                               NA
19:       BD   uvasc                                               NA
20:      MDD   uvasc                                               NA
21:     SCZD   uvasc                                               NA
22:       BD   oligo                                               NA
23:      MDD   oligo                                               NA
24:     SCZD   oligo                                               NA
    disorder context               hit1               hit2        PP3
      <char>  <char>             <char>             <char>      <num>
             PP4
           <num>
 1: 9.212241e-01
 2: 2.425128e-02
 3: 7.604426e-03
 4: 1.373010e-02
 5: 9.824717e-01
 6: 5.597443e-05
 7: 8.244220e-01
 8: 1.486634e-03
 9: 1.137308e-01
10: 8.328080e-01
11: 2.792748e-02
12: 1.251557e-01
13: 8.711159e-01
14: 3.976051e-02
15: 1.006489e-01
16:           NA
17:           NA
18:           NA
19:           NA
20:           NA
21:           NA
22:           NA
23:           NA
24:           NA
             PP4
           <num>
