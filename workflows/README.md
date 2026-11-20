# Breadcrumbs

Bury `.nextflow.log` files in subdir:
```
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/create-spe.log"
```

To specify different .log names for different runs use
```
nextflow -log .nextflow/logs/01_create-spe.log \
run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml
```

Run nextflow from the project/root directory. If running from JHPCE
make sure you're on a node where `sbatch` is available. In this example we are
running the workflow with `-profile test` to use only 3 donors for
testing. `nextflow-test.yaml` supplies the Docker image name and the Apptainer
SIF path. The runtime profile selects which image to use.
```
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml \
-with-report workflows/01_report.html #if want to create report
# if running a continuance on pre-fetched GSE or pre-built spe
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -resume
```

Optionally run and generate resource use report
```
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -with-report workflows/01_report.html
```

Clean up wd and remove empty work subdir
```
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
- Successful Docker run with container, GTF input variable, and curl GEO archive: `spontaneous_avogadro` 
```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test \
  -params-file nextflow.yaml \
  -with-report workflows/reports/Stage01_docker-report.html
#  --geo_archive_resolve_ip 130.14.250.12
```

- Successful Apptainer test run on JHPCE with container, GTF input variable, and curl GEO archive: `boring_picasso`
```bash
srun --mem 1G --pty bash # need to be in compute or transfer node to load nextflow module
module load nextflow/26.10.6
nextflow run workflows/Stage01_create-spe.nf -profile apptainer,test -params-file nextflow-test.yaml \
  -o apptainer-test \
  -with-report workflows/01_apptainer-report.html
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

- Successful Docker run processes 1-3: `fervent_escher`
```bash
nextflow run workflows/Stage02_QC-annotation.nf -profile docker,test \
  -params-file nextflow.yaml 
```


# TO-DO

## Stage 2 QC and annotation

The first stage 2 process, `MAD_OUTLIER_DETECTION`, computes sample- and
slide-level lower 3-MAD flags and detects spatial edges and problem areas.
From the repository root, supply the stage 1 directory containing the paired
`spe_n<n_samples>_assays.h5` and `spe_n<n_samples>_se.rds` files:



For a full run, omit the `test` profile and pass the full stage 1 output
directory. Use `apptainer` in place of `docker` on JHPCE.
Process 1 writes an intermediate metadata table containing the original spot
metadata and the new outlier, edge, and problem-area fields. This table feeds
the review and SpotSweeper processes. Only the expanded SpotSweeper table is
published as `processed-data/Stage02_QC-annotation/colData_QC.csv`.

Process 2, `REVIEW_EDGES_OUTLIER_AREAS`, consumes that table and the stage 1
SPE directory. It publishes `check_edges-outlier-areas.pdf`, with paired
UMI and QC-category plots for each sample flagged for edge artifacts or low
UMI areas. Images are resolved from the stage 1 GEO bundle when available.
Missing images trigger a warning and spots-only plots. If no samples are
flagged, the PDF records that result.

Process 3, `SPOTSWEEPER`, consumes process 1's table and the stage 1 SPE.
It excludes edges, flagged outlier areas, spots with at most 100 UMIs, and
off-tissue spots before computing local outliers for UMI counts, detected
genes, and mitochondrial fractions. It publishes the expanded table to
`processed-data/Stage02_QC-annotation/colData_QC.csv`. Local-outlier values
are NA for spots excluded before SpotSweeper. The script retains the supplied
pre-SpotSweeper `remove_spots` flag.

Optionally set `ignore_edge_samples` to a YAML list of sample IDs to retain
their edge-labeled spots (other exclusion rules still apply):

```yaml
ignore_edge_samples: ['V13F27-338_A1']
```

Omit this key to apply the edge labels to every sample.

**01_create-spe**
- update GEO series matrix to include all donor metadata (age, RIN, seq batch/ round)
