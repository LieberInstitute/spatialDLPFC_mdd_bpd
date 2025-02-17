# Re-run nnSVG on failed samples
After the first run of `nnSVG` using an array job, there were 2 samples that produced the mysterious BiocParallel error. See `processed-data/04_feature_selection/nnSVG_initial-completion_table.txt` for info on the first nnSVG run and the errored samples.
```
Stop worker failed with the error: wrong args for environment subassignment
Error: BiocParallel errors
  0 remote errors, element index:
  6062 unevaluated and other errors
  first remote error:
Execution halted
```
This error previously appeared in more samples and for all other samples it was resolved by removing any gene present in <100 spots (see [nnSVG documentation](https://bioconductor.org/packages/3.21/bioc/vignettes/nnSVG/inst/doc/nnSVG.html#5_Troubleshooting)). After thoroughly evaluating the genes in these two samples, I am confident that the issue is NOT the sparsity of genes.

Running `nnSVG` with `SerialParam()`  based on [this recommendation](https://github.com/sneumann/xcms/issues/627#issuecomment-2084735343) produces a more informative error message, since the genes are not processed simultaneously: (see `nnSVG_min-100_BiocParallel-error-slide333-serial_13819212.log`)
```
Error: BiocParallel errors
  1 remote errors, element index: 1947
  4135 unevaluated and other errors
  first remote error:
Error in BRISC_estimation(coords = coords, y = y_i, x = X, cov.model = "exponential", : c++ error: dpotrf failed
```
Searching for this error produced the following suggestions (all related to spatial coordinates):
- [CRAN update description](https://cran.r-project.org/web/packages/spOccupancy/news/news.html) that stated *"Added in a new check in all spatially-explicit models to see if all the spatial coordinates in the data$coords object were unique, as this is a requirement for spOccupancy spatially-explicit models. In previous versions, this resulted in an error of c++ error: dpotrf failed, or something along those lines, which was a common source of confusion."*
  - I checked and there are no non-unique spatial coordinates or repeated spatialCoord rownames
- [Another suggestion](https://stackoverflow.com/questions/39401698/spbayes-splm-function-with-duplicate-coordinates) related to spatialCoords: *"dup <- duplicated(bef.sp$coords); bef.sp$coords[dup] <- bef.sp$coords[dup] + 1e-3."*
- Another [example](https://github.com/biodiverse/spOccupancy/issues/12) of spatial coordinate repeats but one reference to other potential data formatting issues driving it
- Google AI suggestions: *"This issue typically arises when the input matrix to the dpotrf function (likely a covariance or correlation matrix) is not positive definite. A matrix must be positive definite for the Cholesky decomposition to be computed."*
  - Floating-point precision limitations can lead to numerical instability, especially when dealing with ill-conditioned matrices. If possible, increase the precision of the calculations (e.g., using double instead of float).
    - This is not the case, all rows in `logcounts` are `double` (`table(apply(logcounts(tmp2), 1, typeof))`)
  - Problems in the input data, such as highly correlated variables or outliers, can result in a non-positive definite matrix. 

Both these samples (`V13B23-302_B1` and `V13B23-333_C1`) previously ran with when subsetting per-slide spe and using `conda_R/devel` module. Attempted to replicate the archived code (`03-tmp_nnSVG-slide333_replicate-per-sample`) while using `conda_R/4.4.x` and found that package versions were not back compatible (see log `nnSVG_slide333_per-sample-replicate_13868217.log`). Looked at version differences for necessary packages and found enough differences to try running these samples with `conda_R/devel`.
```
$ module load conda_R/4.4.x
> packageVersion("BiocParallel")
[1] '1.40.0'
> packageVersion("nnSVG")
[1] ‘1.10.2
> packageVersion("HDF5Array")
[1] ‘1.34.0’
> packageVersion("DelayedArray")
[1] ‘0.32.0’
> packageVersion("BRISC")
[1] ‘1.0.6’
> packageVersion("SpatialExperiment")
[1] ‘1.16.0’

$ module load conda_R/devel
> packageVersion("BiocParallel")
[1] ‘1.41.0’
> packageVersion("nnSVG")
[1] ‘1.11.0’
> packageVersion("HDF5Array")
[1] ‘1.35.13’
> packageVersion("DelayedArray")
[1] ‘0.33.6’
> packageVersion("BRISC")
[1] ‘1.0.6’
> packageVersion("SpatialExperiment")
[1] ‘1.17.0’
```

First try running the same nnSVG script but now in the `conda_R/devel` environment (`03-dev_nnSVG_per-sample_RERUN.sh`) produced this error for both jobs:
(see logs `standard-nnSVG_per-sample_min-100_RERUN-devel_13882012_1.log` and `standard-nnSVG_per-sample_min-100_RERUN-devel_13880316_4.log`)
```
Error in reducer$value.cache[[as.character(idx)]] <- values :
  wrong args for environment subassignment
Calls: nnSVG ... .bploop_impl -> .collect_result -> .reducer_add -> .reducer_add
```

Currently trying to run again but forcing the job to occur on only 1 node. Next step would be to build single-sample spe in the dev environment and try again. Another try could be to update select packages in `conda_R/4.4.x` module (start with BiocParallel to dev version?)

### Log output
Looking at the logs is really clunky because of the verbose output so I use the following code in `R` to navigate/pull information:
```
x <- readLines("code/04_feature_selection/logs/nnSVG_min-100_BiocParallel-error-slide333-serial_13819212.log")
x[(length(x)-10):length(x)] #takes you straight to the end/error
x[19:25] #gives you the sample number and the next couple lines give the # of genes filtered and the spe dims
x[grep("333_",x)] #If it has completed the sample number will be in save lines, if not the last line with the sample will read that the model started
length(grep("BRISC model fit with 4951 observations",x)) #gives you the number of genes that ran before it errored out
```
# spoon/ weighted_nnSVG
Looked into using spoon/weighted_nnSVG to further adjust for the mean-rank relationship but found that the difference between weighted_nnSVG and standard nnSVG wasn't obvious. Since we still would've wanted to select top SVGs from expression quantiles manually, decided to stick with standard nnSVG. In the `archive/` dir of `code/04_feature_selection/` and `processed-data/04_feature_selection/` there are many files related to running spoon/weighted_nnSVG.
