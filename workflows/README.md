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
- pending Apptainer run
```bash
srun --mem 1G --pty bash # need to be in compute or transfer node to load nextflow module
module load nextflow # current v24.10.5 has lots of syntax differences from 26.10.0 so asking for update
nextflow run workflows/01_create-spe.nf -profile apptainer,test -params-file nextflow-test.yaml \
  -o nextflow_outputs/apptainer-test 
```

# TO-DO

**01_create-spe**
- write options for reference genome file: currently pulling from features.csv which
option I want to keep but also want to be able to set an absolute path to a desired
gtf file. Final product will be an input param stored in .yaml clarifying the 
reference genome input
