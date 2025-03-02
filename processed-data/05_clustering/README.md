# Description of different SVG filters used

## n=1053
- padj<.05 and max rank across all samples <= 500
- padj<.05 in at least 10 samples
- padj<.05 in at least 33% of samples
- Not a batch effect gene (n=43)
- Most highly expressed (decile==10) qualifying genes removed with filter of round(spcov,2)>=0.4
## n=1050
- padj<.05 and max rank across all samples <= 500
- padj<.05 in at least 10 samples
- padj<.05 in at least 33% of samples
- Not a batch effect gene (n=43)
- Most highly expressed (decile==10) qualifying genes removed with filter of round(spcov,2)>=0.4
- Remove the 3 remaining RPS | RPL genes
## n=1051
When I followed the guidance of re-calculating the adj. p value from nnSVG *after* filtering out genes present in <500 spots, there was 1 additional gene in the final SVG list: POLR2F. Unfortunately the inclusion of this single gene also changes the PRECAST clusters, very much like the inclusion of the 3 RPS|RPL genes.

COME BACK AND EVALUATE THIS AFTER MAKING HPC FIGURES. Since I'm going to try BayesSpace anyway, maybe it won't matter with that method??
