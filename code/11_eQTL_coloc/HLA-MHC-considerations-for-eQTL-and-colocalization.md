# HLA/MHC Considerations for eQTL and Colocalization

## Local Evidence

The current eQTL workflow maps genome-wide cis-eQTLs with tensorQTL and compares significant eGenes with the broad PRECAST plus Seurat DEG union. HLA/MHC signals are present, but the DEG x eGene intersection is concentrated in a small number of repeated HLA genes.

Current local counts:

- Broad DEG union: 817 genes.
- HLA genes in the broad DEG union: 4 (`HLA-B`, `HLA-DPA1`, `HLA-DRA`, `HLA-DRB1`).
- Seurat significant lead eGenes, `qval < 0.05`: 1435 unique genes.
- Seurat eGene-DEG overlap: 79 unique genes.
- Seurat HLA eGenes: 9 (`HLA-A`, `HLA-B`, `HLA-C`, `HLA-DMA`, `HLA-DMB`, `HLA-DPB1`, `HLA-DQB1`, `HLA-DRB1`, `HLA-DRB5`).
- Seurat HLA eGene-DEG overlap: 2 genes (`HLA-B`, `HLA-DRB1`).
- Seurat MHC-region eGenes using hg38 `chr6:25000000-35000000`: 21 unique genes.
- Custom-cluster results are similar: 1239 unique eGenes, 73 DEG-overlap eGenes, 9 HLA eGenes, 2 HLA DEG-overlap eGenes, and 20 MHC-region eGenes.

The apparent prominence of HLA genes in summary tables is driven partly by recurrence across contexts and splits. In the DEG x eGene intersection, the unique HLA contribution is mainly `HLA-B` and `HLA-DRB1`.

## Literature Screening

The HLA/MHC region is a special case for genetic mapping. D'Antonio et al. describe the MHC as highly polymorphic, gene dense, and marked by strong linkage disequilibrium. Their MHC-focused eQTL analysis found that HLA-type haplotypes provided greater power than single-variant eQTL analysis for MHC-region eGenes. Source: D'Antonio et al., eLife 2019, "Systematic genetic analysis of the MHC region reveals mechanistic underpinnings of HLA type associations with disease": https://elifesciences.org/articles/48476v1

HLA expression quantification from standard RNA-seq references can be biased. Aguiar et al. note that short reads from highly polymorphic HLA genes can fail to map to the reference genome or map ambiguously across paralogues. Their personalized HLA pipeline was designed to reduce those biases and improve HLA expression and eQTL mapping. Source: Aguiar et al., PLOS Genetics 2019, "Expression estimation and eQTL mapping for HLA genes with a personalized pipeline": https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008091

Single-cell HLA eQTL work also treats HLA as technically distinct. Kang et al. used personalized reference genomes to mitigate technical confounding and identified cell-type-specific cis-eQTLs for classical HLA genes. Source: Kang et al., Nature Genetics 2023, "Mapping the dynamic genetic regulatory architecture of HLA genes at single-cell resolution": https://www.nature.com/articles/s41588-023-01586-6

Colocalization in complex LD regions requires caution. Giambartolomei et al. emphasize that observing a shared associated SNP is not sufficient to establish colocalization, because distinct causal variants can be in linkage disequilibrium. Source: Giambartolomei et al., Human Molecular Genetics 2015, "Integration of disease association and eQTL data using a Bayesian colocalisation approach highlights six candidate causal genes in immune-mediated diseases": https://academic.oup.com/hmg/article/24/12/3305/621728

Some modern colocalization workflows exclude the HLA/MHC region from main analyses because of LD complexity. For example, Soskic et al. excluded hg19/GRCh37 `chr6:25000000-35000000` in a cross-disorder immune-disease colocalization workflow. Source: Soskic et al., Nature Communications 2023, "Cross-disorder genetic analysis of immune diseases reveals distinct gene associations that converge on common pathways": https://www.nature.com/articles/s41467-023-38389-6

Broader eQTL resources also highlight challenges relevant to colocalization, including multiple independent signals and the need for accurate LD information. Source: Kerimov et al., Nature Genetics 2021, "A compendium of uniformly processed human gene expression and splicing quantitative trait loci": https://www.nature.com/articles/s41588-021-00924-w

## Workflow Implications

The presence of HLA eGenes does not by itself indicate that the tensorQTL cis-eQTL workflow is invalid. HLA eQTLs are biologically plausible and expected, and the broader MHC-region eGene count is not globally inflated relative to all tested genes in this dataset.

The main workflow should remain genome-wide cis-eQTL discovery. HLA/MHC handling is best treated as a reporting and interpretation sensitivity layer unless HLA biology becomes a central claim.

The recommended reporting approach is:

- Keep tensorQTL mapping unchanged.
- Add HLA/MHC flags to row-level eQTL tables.
- Report DEG x eGene and coloc-facing summaries with all rows, excluding HLA genes, and excluding MHC-region genes or lead variants.
- Treat HLA/MHC colocalization results as exploratory unless a dedicated HLA-aware analysis is performed.
- Avoid strong causal language for HLA/MHC colocalization from standard `coloc.abf`, especially in the presence of psychiatric GWAS signals and dense LD.

## Recommended Language

Suggested neutral wording:

"The HLA/MHC region was retained in the primary cis-eQTL discovery workflow. Because the region is highly polymorphic, gene dense, and characterized by complex linkage disequilibrium, HLA/MHC results are additionally reported with sensitivity filters excluding HLA genes and excluding the broader MHC interval. Colocalization results involving HLA/MHC loci should be interpreted as exploratory unless supported by HLA-aware quantification, HLA typing, or fine-mapping designed for this region."

Specialized HLA typing or personalized HLA expression quantification is outside the current implementation scope. Such analyses would be appropriate only if HLA/MHC regulation becomes a primary biological claim.
