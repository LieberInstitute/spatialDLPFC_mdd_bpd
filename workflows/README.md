# Breadcrumbs

Bury `.nextflow.log` files in subdir:
```bash
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/create-spe.log"
```

To specify different .log names for different runs use
```bash
nextflow -log .nextflow/logs/01_create-spe.log \
run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml
```

Run nextflow from the project/root directory. If running from JHPCE
make sure you're on a node where `sbatch` is available. In this example we are
running the workflow with `-profile test` to use only 3 donors for
testing. `nextflow-test.yaml` supplies the Docker image name and the Apptainer
SIF path. The runtime profile selects which image to use.
```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml \
-with-report workflows/01_report.html #if want to create report
# if running a continuance on pre-fetched GSE or pre-built spe
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -resume
```

Optionally run and generate resource use report
```bash
nextflow run workflows/Stage01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -with-report workflows/01_report.html
```

Clean up wd and remove empty work subdir
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

- Apptainer full run on JHPCE: ``
```bash
srun --mem 4G --pty bash # previously had nextflow success but OOM runner with 1GB
module load nextflow/26.10.6
export NXF_JVM_ARGS="-Xms256m -Xmx2g"
nextflow run workflows/Stage02_QC-annotation.nf -profile apptainer -params-file nextflow.yaml \
  -with-report workflows/reports/Stage02_apptainer-full-report.html \
  --spe_dir nextflow_outputs/01_create-spe
```

# TO-DO

**Stage01_create-spe**
- update GEO series matrix to include all donor metadata (age, RIN, seq batch/ round)

**Stage02_QC-annotation**
- make QC specific readme: explain ignore_edge_samples and n=120 vs n=119 differences?