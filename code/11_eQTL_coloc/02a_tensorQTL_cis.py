#!/usr/bin/env python
# ruff: noqa: E402
## this should use the activated environment where tensorQTL and rpy2 are installed

## hack to fix libstdc++ conflict with R
import rpy2.robjects.packages as rp
rp.importr("qvalue")
## - now import the rest of packages that might pull the system libstdc++ version
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis
#import time
import os
import os.path as op
import sys
import glob
import re
import argparse
from pathlib import Path
print(f'PyTorch {torch.__version__}')
cuda_ok = torch.cuda.is_available()
cuda_dev = torch.cuda.get_device_name(torch.cuda.current_device()) if cuda_ok else "CPU"
print(f'CUDA available: {cuda_ok} ({cuda_dev})')
print(f'Pandas {pd.__version__}')
print(f'tensorqtl {tensorqtl.__version__}')

repo_root = Path(__file__).resolve().parents[2]
default_in_dir = str(repo_root / 'processed-data' / '11_eQTL_coloc' / 'tqtl_in')
default_out_dir = str(repo_root / 'processed-data' / '11_eQTL_coloc' / 'tqtl_out')
default_plink = str(repo_root / 'processed-data' / '00_genotypes' / 'plink2' / 'merged_maf05')

parser = argparse.ArgumentParser(
    description='Run tensorQTL cis-mapping for one prepared dataset ID'
)
parser.add_argument('dsname', help='Dataset prefix matching <input-dir>/<dsname>.gene.expr.bed.gz')
parser.add_argument('--input-dir', default=default_in_dir, help='Directory containing *.gene.expr.bed.gz and *.gene.covars.txt')
parser.add_argument('--output-dir', default=default_out_dir, help='Directory for tensorQTL outputs')
parser.add_argument('--plink-prefix', default=default_plink, help='PLINK2 prefix for genotype data (.pgen/.psam/.pvar)')
cli = parser.parse_args()

outdir = cli.output_dir
os.makedirs(outdir, exist_ok=True)

in_dir = cli.input_dir
dsname = cli.dsname

run_nominal_all = True  # set False to only run nominal for significant eGenes
egenes_qval = 0.10

plink = cli.plink_prefix
#plink='genotypes/mdd_maf01'
mapping_file=None ## comment the line below if prep_tensorQTL already used geno2rna mapping
#mapping_file='data/geno2rna.tab' # genotype IDs will be changed to their RNAseq mappings in the 2nd column
#mapping_file='genotypes/mdd_geno2rna.tab'
## if not None, mapping_file must be the same with the one used for 02_prep_tensorQTL_RSE_input.R
col_map = None
if mapping_file is not None and op.exists(mapping_file):
    print(f"mapping genotype IDs to RNAseq sample IDs based on mapping file: {mapping_file}")
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["gt_id", "r_id"])
    col_map = dict(zip(mapping_df['gt_id'], mapping_df['r_id']))
if op.exists(plink+'.pgen'): # checks for plink2 pgen file for genotype data
    print("Plink2 pgen file found.")
else:
    raise Exception("Plink2 "+plink+".pgen not found!")

exprfiles=glob.glob(in_dir+'/'+dsname+'.*.expr.bed.gz')
features = [os.path.basename(file).split('.')[1] for file in exprfiles]
exprdict = {os.path.basename(file).split('.')[1]: file for file in exprfiles}
for feat in features:
    expres = exprdict[feat]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    if not op.exists(covar):
       raise FileNotFoundError(f"Covars file {covar} not found!")
print("Features found: ",features)
#print(" Change the `features` array below to set the features to be processed:")
## CHANGE here and uncomment, if needed, in order to select the features to process:
#features = ['gene']


## loading genotype data
#pr = genotypeio.PlinkReader(plink)
## use plink2 pgen format
pgr = pgen.PgenReader(plink)
genotype_df = pgr.load_genotypes()
#variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
variant_df = pgr.variant_df
print("Genotype dimensions: ", end='')
print(genotype_df.shape)
# if a mapping file is given
# Check if all genotype_df column names are in the mapping file
#all_columns_found = genotype_df.columns.isin(mapping_df["current_column"]).all()

if col_map is not None:
  # Rename columns using the dictionary
  genotype_df.rename(columns=col_map, inplace=True)
## Now genotype_df contains the updated column names
#print(genotype_df.iloc[:5, :7])  # Display the first few rows

## run tensorQTL for each feature
def fixBrNums(col):
    return re.sub(r'^Br(\d\d\d)$', r'Br0\1', col)

## fix Brnums with 3 digits, if any
genotype_df.columns = [fixBrNums(col) for col in genotype_df.columns]

##  etc. Fix chromosomes (add the "chr" prefix) if needed:
variant_df.chrom = variant_df.chrom.astype(str)
if not variant_df.chrom.str.startswith('chr').all():
    variant_df.chrom = variant_df.chrom.map(
        lambda chrom: chrom if chrom.startswith('chr') else 'chr' + chrom
    )
## select chromosomes - to make sure we have the same chromosomes in our data for each expression dataset
variant_chrom = set(variant_df.chrom)

for feat in features:
    print(f" Processing {feat} features for dataset {dsname} ..")
    tag = dsname+'.'+feat
    expres = exprdict[feat]
    covar=expres.replace("expr.bed.gz", "covars.txt")
    covariates_df = pd.read_csv(covar, sep='\t', index_col=0).T
    ## fix Brnums with 3 digits, just in case
    covariates_df.index = [fixBrNums(idx) for idx in covariates_df.index]

    numcovars = covariates_df.shape[1]
    print(f"Using {numcovars} covariates: ", covariates_df.columns.tolist())

    #print(covariates_df.iloc[:5, :5])
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expres)
    print("Phenotype dimensions:", end='')
    print(phenotype_df.shape)
    ## fix Brnums with 3 digits
    phenotype_df.columns = [fixBrNums(col) for col in phenotype_df.columns]

    ## use the same chromosome set
    express_chrom = set(phenotype_pos_df.chr)
    assert len(variant_chrom.intersection(express_chrom))>0
    if express_chrom - variant_chrom:
        chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
        if chrom_filter.sum()<phenotype_df.shape[0]:
          phenotype_df = phenotype_df[chrom_filter]
          phenotype_pos_df = phenotype_pos_df[chrom_filter]
    ## make sure we keep only the genotypes for the same expression samples
    cols=phenotype_df.columns.tolist()
    gcols=genotype_df.columns.tolist()
    if phenotype_df.columns.duplicated().any():
        raise RuntimeError("Duplicate sample IDs detected in expression data.")
    ## which genotypes are in the expression data? print them here if any are missing and stop the process
    missing_geno = set(cols) - set(gcols)
    if missing_geno:
        print(f"Genotypes missing for {len(missing_geno)} samples: {missing_geno}")
        raise RuntimeError("Missing genotypes detected. Aborting...")
    ## show if any genotypes are not in the expression data (non-fatal):
    missing_expr = set(gcols) - set(cols)
    if missing_expr:
        print(f"These genotypes have no expression data given: {missing_expr}")
    geno_df=genotype_df.loc[:, cols]
    covar_index = covariates_df.index.tolist()
    if covariates_df.index.duplicated().any():
        raise RuntimeError("Duplicate sample IDs detected in covariates.")
    missing_covar = set(cols) - set(covar_index)
    if missing_covar:
        print(f"Covariates missing for {len(missing_covar)} samples: {missing_covar}")
        raise RuntimeError("Missing covariates detected. Aborting...")
    extra_covar = set(covar_index) - set(cols)
    if extra_covar:
        print(f"These covariate rows have no expression data given: {extra_covar}")
    covariates_df = covariates_df.loc[cols]
    non_numeric = [
        col for col in covariates_df.columns
        if not pd.api.types.is_numeric_dtype(covariates_df[col])
    ]
    if non_numeric:
        raise RuntimeError(f"Non-numeric covariates detected: {non_numeric}")
    if covariates_df.isna().any().any():
        raise RuntimeError("Missing values detected in covariates.")
    ## run tensorQTL, map_cis:
    #cis.map_nominal(geno_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
    #            maf_threshold=0.05, window=500000, output_dir= outdir, verbose=False)
    ## run permutation-based mapping with map_cis (takes longer)
    ## cis.map_cis is LD aware and it only outputs the most significant eQTL per gene
    cis_out = os.path.join(outdir, f"{tag}.map_cis.tab.gz")
    if os.path.exists(cis_out):
        cis_df = pd.read_csv(cis_out, sep='\t', index_col=0)
    else:
        print(f" .. running map_cis for {feat} ..")
        cis_df = cis.map_cis(geno_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
            group_s=None, maf_threshold=0.05, beta_approx=True, nperm=10000,
            window=1000000, random_tiebreak=False, seed=118, verbose=False
        )
        tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85, fdr=0.05) ## this adds qval columns to cis_df
        cis_df.to_csv(cis_out, sep='\t', compression='gzip', index=True)
    ## run map_nominal
    if run_nominal_all:
        print(f" .. running map_nominal for {feat}, all phenotypes ..")
        ph_sub = phenotype_df
        phpos_sub = phenotype_pos_df
        cis.map_nominal(
            geno_df, variant_df, ph_sub, phpos_sub, prefix=tag,
            covariates_df=covariates_df,
            maf_threshold=0.05, window=1000000, run_eigenmt=True,
            output_dir=outdir, write_top=False, write_stats=True, verbose=False
        )
    else:
        egenes = cis_df.index[cis_df['qval'] < egenes_qval]
        if len(egenes) > 0:
            print(f" .. running map_nominal for {feat}, {len(egenes)} eGenes ..")
            ph_sub = phenotype_df.loc[egenes]
            phpos_sub = phenotype_pos_df.loc[egenes]
            cis.map_nominal(
                geno_df, variant_df, ph_sub, phpos_sub, prefix=tag,
                covariates_df=covariates_df,
                maf_threshold=0.05, window=1000000, run_eigenmt=True,
                output_dir=outdir, write_top=False, write_stats=True, verbose=False
            )
        else:
            print(f"No phenotypes with q<{egenes_qval:.2f}; skipping map_nominal.")
    indep_out = os.path.join(outdir, f"{tag}.map_independent.txt.gz")
    ## run map_independent (skip safely if no significant phenotypes at selected FDR)
    indep_fdr = 0.05
    if 'qval' not in cis_df.columns:
        print(f"Skipping map_independent for {feat}: 'qval' column not found in cis_df.")
    else:
        n_sig = int((cis_df['qval'] <= indep_fdr).sum())
        if n_sig == 0:
            print(f"Skipping map_independent for {feat}: no significant phenotypes at FDR <= {indep_fdr}.")
        else:
            print(f" .. running map_independent for {feat} ({n_sig} significant phenotypes) ..")
            try:
                indep_df = cis.map_independent(
                    geno_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df,
                    fdr=indep_fdr, verbose=False
                )
                indep_df.to_csv(indep_out, sep='\t', compression='gzip', index=True)
            except ValueError as e:
                if "No significant phenotypes at FDR" in str(e):
                    print(f"Skipping map_independent for {feat}: {e}")
                else:
                    raise
    ## Dx interaction scan, top only
    # Support both legacy single-column DX and one-hot columns from model.matrix (e.g. DxMDD, DxMDD_OUD).
    dx_candidates = [c for c in ['DXSCZ', 'DX', 'Dx', 'dx'] if c in covariates_df.columns]
    dx_prefix_cols = [c for c in covariates_df.columns if re.match(r'(?i)^dx', c)]
    for c in dx_prefix_cols:
        if c not in dx_candidates:
            dx_candidates.append(c)

    valid_dx_cols = []
    for c in dx_candidates:
        x = pd.to_numeric(covariates_df[c], errors='coerce')
        if x.isna().any():
            print(f"Skipping interaction covariate {c}: non-numeric values detected.")
            continue
        ux = set(x.unique())
        if ux.issubset({0, 1}):
            valid_dx_cols.append(c)
        else:
            print(f"Skipping interaction covariate {c}: values are not strictly 0/1 -> {ux}")

    if len(valid_dx_cols) == 0:
        print("No usable Dx interaction covariate found; skipping interaction scan.")
    else:
        for dx_col in valid_dx_cols:
            print(f" .. running map_nominal interaction for {feat} ({dx_col}) ..")
            covar_ix = covariates_df.drop(columns=[dx_col])
            dx = pd.to_numeric(covariates_df[dx_col], errors='coerce')
            cis.map_nominal(
                geno_df, variant_df, phenotype_df, phenotype_pos_df, prefix=tag + f".{dx_col}INT",
                covariates_df=covar_ix,
                interaction_df=dx.to_frame('Dx'), maf_threshold_interaction=0.05,
                run_eigenmt=True, group_s=None, window=1000000,
                output_dir=outdir, write_top=True, write_stats=False, verbose=False
            )
print("All done.")
