# Breadcrumbs

Bury `.nextflow.log` files in subdir:
```
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/create-spe.log"
```

To specify different .log names for different runs use
```
nextflow -log .nextflow/logs/01_create-spe.log \
run workflows/01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml
```

Run nextflow from the project/root directory. If running from JHPCE
make sure you're on a node where `sbatch` is available. In this example we are
running the workflow with `-profile test` to use only 3 donors for
testing. `nextflow-test.yaml` supplies the Docker image name and the Apptainer
SIF path. The runtime profile selects which image to use.
```
nextflow run workflows/01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml \
-with-report workflows/01_report.html #if want to create report
# if running a continuance on pre-fetched GSE or pre-built spe
nextflow run workflows/01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -resume
```

Optionally run and generate resource use report
```
nextflow run workflows/01_create-spe.nf -profile docker,test -params-file nextflow-test.yaml -with-report workflows/01_report.html
```

Clean up wd and remove empty work subdir
```
nextflow clean -before <run_name> -f && find work -type d -empty -delete
```

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
- Successful Docker run with container, GTF input variable, and curl GEO archive: `mighty_bhaskara` 
```
nextflow run workflows/01_create-spe.nf -profile docker,test \
  -params-file nextflow-test.yaml \
  --geo_archive_resolve_ip 130.14.250.12
```

- Successful Apptainer run on JHPCE with container, GTF input variable, and curl GEO archive: `boring_picasso`
```bash
srun --mem 1G --pty bash # need to be in compute or transfer node to load nextflow module
module load nextflow/26.10.6
nextflow run workflows/01_create-spe.nf -profile apptainer,test -params-file nextflow-test.yaml \
  -o nextflow_outputs/apptainer-test \
  -with-report workflows/01_apptainer-report.html
```


# TO-DO

**01_create-spe**
- update GEO series matrix to include all donor metadata (age, RIN, seq batch/ round)
