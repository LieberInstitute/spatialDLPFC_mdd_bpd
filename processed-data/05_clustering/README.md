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
