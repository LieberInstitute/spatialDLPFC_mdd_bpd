# Breadcrumbs

Bury `.nextflow.log` files in subdir:
```
mkdir -p .nextflow/logs
export NXF_LOG_FILE="$PWD/.nextflow/logs/pseudobulk.log"
```

Clean up wd and remove empty work subdir
```
nextflow clean -before <run_name> -f && find work -type d -empty -delete
```

Run nextflow
```
nextflow run main.nf -params-file nextflow-test.yaml
```
