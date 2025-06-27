### Associated with code/05_clustering/Seurat
As an alternative to unsupervised clustering, we explored using `Seurat` to transfer labels from snRNA-seq DLPFC data to our Visium SRT spots.

After consulting with Kinnary (who explored this for the dACC project), we looked at using published PsychENCODE datasets.

Following her pipeline, we accessed the `V5` data downloaded by Nick from here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/`

The `SZBDMulti-Seq_annotated.h5ad` file was loaded directly from that location and processed as described in `code/05_clustering/Seurat/01_adata_processing.ipynb`. As a backup, I copied the source H5AD file to this location.
