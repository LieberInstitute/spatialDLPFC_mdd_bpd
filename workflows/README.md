# Breadcrumbs

Bury `.nextflow.log` files in subdir:
```
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/create-spe.log"
```

To specify different .log names for different runs use
```
nextflow -log .nextflow/logs/01_create-spe.log \
run workflows/01_create-spe.nf -profile test
```

Run nextflow from the project/root directory. In this example we are
running the workflow with `-profile test` to use only 3 donors for
testing.
```
nextflow run workflows/01_create-spe.nf -profile test
# if running a continuance on pre-fetched GSE or pre-built spe
nextflow run workflows/01_create-spe.nf -profile test -resume
```

Optionally run and generate resource use report
```
nextflow run workflows/01_create-spe.nf -profile test -with-report workflows/01_report.html
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

# TO-DO

**01_create-spe**
- write options for reference genome file: currently pulling from features.csv which
option I want to keep but also want to be able to set an absolute path to a desired
gtf file
