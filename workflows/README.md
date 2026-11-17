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

Run nextflow from the project/root directory. In this example we are
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

### Apptainer test on JHPCE

Run from the repository root on a JHPCE node where `sbatch` is available.
Nextflow submits each task as a Slurm job. Load the existing Apptainer module
and use the `apptainer,test` profiles together:

```bash
module load apptainer/1.5.3
nextflow run workflows/01_create-spe.nf -profile apptainer,test -params-file nextflow-test.yaml --timeout 600
```

The profile uses the `mbv-r-jhpce.sif` from build `20261112T171336Z` in the
original JHPCE checkout. If the SIF lives elsewhere, set `r_apptainer_image`
in the params file to its absolute path. Test outputs go to
`nextflow_outputs/test/01_create-spe/`. The work filesystem needs room for
the 4.4 GB GEO archive and its extracted contents.

Clean up wd and remove empty work subdir
```
nextflow clean -before <run_name> -f && find work -type d -empty -delete
```

# Tests completed

**01_create-spe**
- successful Docker run: `deadly_waddington` 
```
nextflow run workflows/01_create-spe.nf -profile test -with-report workflows/01_report.html --timeout 600
```
- successful Docker run (with container revisions): `irreverent_celsius`
```
nextflow run workflows/01_create-spe.nf -profile docker,test \
  -params-file nextflow-test.yaml
```

# TO-DO

**01_create-spe**
- write options for reference genome file: currently pulling from features.csv which
option I want to keep but also want to be able to set an absolute path to a desired
gtf file. Final product will be an input param stored in .yaml clarifying the 
reference genome input
