# Breadcrumbs
All commands executed in the cloned project/root directory. If running from JHPCE
make sure you're on a node where `sbatch` is available.

Bury `.nextflow.log` files in subdir:
```bash
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/.nextflow.log"
```

Pick a nextflow entry point to run a portion of the MBv analysis. There are 4 entry
points for this project, separating the analysis into 4 stages. 

Option to build SPE from FASTQ (Stage 1)
Alternatively can construct unfiltered SPE from GSE333974 supplemental files.
Either of those SPE can be input to reproduce QC and annotation with Stage 2 workflow.
Alternatively, use the provided columns in GSE333974 supplemental file (observations).
Provide these to Stage 3 pseudobulk and DE to reproduce statistical model results.

```
MBv analysis pipeline
│
├── Stage 1: Create SPE
│     └── code/02_build_spe
│
├── Stage 2: QC, feature selection & annotation
│     ├── code/03_QC
│     ├── code/04_feature_selection
│     └── code/05_clustering
│
├── Stage 3: Pseudobulk & differential expression
│     ├── code/06_pseudobulk
│     ├── code/07-1_covariate_sensitivity
│     └── code/07_dx_DE
│
└── Stage 4: Gene regulatory networks *
      ├── code/09_DEG_GRN
      └── code/10_SCENIC
```
* Stage 4 is not perfectly reproducible because every run with GRNBOOST2 produces a unique result.

Run a stage of analysis by calling the entry point nextflow file and supplying options
for the following required parameters:
- YAML location: many things, importantly the docker images and apptainer SIF path
- profile: Select docker or apptainer depending on your setup

```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow.yaml
```

There are also some optional parameters, all of which are helpful when troubleshooting and
establishing a workflow.
- resume: pickup from cached checkpoint
- 'test' profile overrides default n=119 samples, running on a subset of 7 samples
- with-report: I use to track usage so that you can anticipate memory and time requirements

```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow.yaml \
-with-report workflows/reports/Stage01_docker-test.html \
-resume
```

Clean up wd (all before specific run) and remove empty work subdir
```bash
nextflow clean -before <run_name> -f && find work -type d -empty -delete
```



### Published output layout

The default output root is the launch directory. The `test` profile uses
`test/` beneath that directory. Each stage publishes analysis files under
`processed-data/<stage_subdir>/` and figures under `plots/<stage_subdir>/`.

For example, the stage 2 test publishes its CSV to
`test/processed-data/Stage02_QC-annotation/` and its PDF to
`test/plots/Stage02_QC-annotation/`. A custom `-o my-run` places both
trees beneath `my-run/`. Supply stage 2's `--spe_dir` from the matching
stage 1 `processed-data/Stage01_create-spe/` directory.

### GEO transfer instability
GEO download speeds and outcomes can vary substantially. To give users 
better flexibility to troubleshoot connections, the large (4GB) tar file
containing spaceranger output for all 119 donors is now retrieved in a
separate workflow step. 

`DOWNLOAD_GEO_ARCHIVE` fetches the tar with curl before `FETCH_GEO` extracts
the selected samples. The archive stays in Nextflow's work directory so
`-resume` can reuse a completed download. Each transfer attempt uses
`--timeout` as its time limit. The download tries up to three times.

If the GEO hostname gives inconsistent responses, you can temporarily add
`--geo_archive_resolve_ip <verified-IP>` to the launch command. Confirm that
the address works from JHPCE before using it. By default, curl uses DNS.

# Tests completed

**01_create-spe**
- Successful Docker run (7 sample) with container, GTF input variable, and curl GEO archive: `condescending_linnaeus` 
```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test \
  -params-file nextflow.yaml \
  -with-report workflows/reports/Stage01_docker-report.html
#  --geo_archive_resolve_ip 130.14.250.12
```

- Successful Apptainer test run with corrected GEO metadata (platform-specific series matrices): `pedantic_austin`
```bash
srun --mem 4G --pty bash # previously had nextflow success but OOM runner with 1GB
module load nextflow/26.10.6
export NXF_JVM_ARGS="-Xms256m -Xmx2g"
export NXF_LOG_FILE="$PWD/.nextflow/logs/.nextflow.log"
nextflow run workflows/Stage01_create-spe.nf -profile apptainer,test -params-file nextflow.yaml \
  -with-report workflows/reports/Stage01_apptainer-test-report.html \
  -w "$MYSCRATCH/MBv-work" -resume
```

- Successful Apptainer full run on JHPCE with container, GTF input variable, and curl GEO archive: `awesome_leibniz`
```bash
srun --mem 4G --pty bash # previously had nextflow success but OOM runner with 1GB
module load nextflow/26.10.6
export NXF_JVM_ARGS="-Xms256m -Xmx2g"
nextflow run workflows/Stage01_create-spe.nf -profile apptainer -params-file nextflow-test.yaml \
  -with-report workflows/01_apptainer-full-report.html \
  -resume
```

**02_QC-annotation**

- Successful Docker run (7 sample) processes 1-5: `grave_kalman`
```bash
nextflow run workflows/Stage02_QC-annotation.nf -profile docker,test \
  -params-file nextflow.yaml 
```

- Apptainer test run on JHPCE with corrected GEO metadata: `condescending_newton` 
```bash
srun --mem 4G --pty bash # previously had nextflow success but OOM runner with 1GB
module load nextflow/26.10.6
export NXF_JVM_ARGS="-Xms256m -Xmx2g"
nextflow run workflows/Stage02_QC-annotation.nf -profile apptainer,test -params-file nextflow.yaml \
  -with-report workflows/reports/Stage02_apptainer-test-report.html \
  --spe_dir test/processed-data/Stage01_create-spe \
  -w "$MYSCRATCH/MBv-work" -resume fa92fb74-8832-4855-8310-975e6a5be892
```

# TO-DO

- container distribution
- load stage-specific resource configurations only for their corresponding entry points to avoid the config warnings

**Stage01_create-spe**
- update GEO series matrix to include all donor metadata (age, RIN, seq batch/ round)

**Stage02_QC-annotation**
- make QC specific readme: explain ignore_edge_samples and n=120 vs n=119 differences?
- QC pdf edge check and 2 sample handling
- nnSVG runs should complete within 2 hours. if much longer times, try again using resume
- make sure different seed is well enough explained