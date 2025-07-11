# Label transfer from snRNA-seq

To better identify L4 and GABA-enriched spots we used Seurat's label transfer method to transfer annotations from DLPFC snRNA-seq data to the MBv SRT dataset.

Kinnary had explored something similar for the dACC project and had compiled several DLPFC snRNA-seq datasets on the server already.

```
https://github.com/LieberInstitute/spatialdACC/blob/main/code/18_PsychENCODE_NMF/00_downloading.py
https://github.com/LieberInstitute/spatialdACC/blob/main/code/18_PsychENCODE_NMF/01_pseudobulking.R
https://github.com/LieberInstitute/spatialdACC/blob/main/code/18_PsychENCODE_NMF/02_preprocessing.R
```
She worked with PsychENCODE datasets, compiled by the [brainSCOPE Resource](https://brainscope.gersteinlab.org/), downloaded by Nick Eagles via the [Synapse data repository](https://www.synapse.org/Synapse:syn51111084.5/datasets/). 

Those downloads are found here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/`

The `SZBDMulti-seq` files used in this analysis are backed up in this project dir: `raw-data/PsychENCODE_backup`

### Selecting DLPFC snRNA-seq data

Scripts exploring DLPFC snRNA-seq datasets are saved locally by JT named `exploring_different_snRNA-seq_datasets.R`.

I first explored using the `spatialDLPFC` n=30 snRNA-seq data set but found the nuclei quality too low.

I then looked at all of the PEv5 datasets by accessing the metadata here: `https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt`. All of these datasets had been made available by Nick. I examined the number of control samples, distribution of control ages, and number of nuclei per control sample. From this I selected the CMC and SZBDMulti-seq datasets. Looking at the library size and number of detected genes showed the SZBDMulti-seq dataset to be of higher and more consistent quality. It has the bonus of also containing high quality bipolar disorder snRNA-seq data.
