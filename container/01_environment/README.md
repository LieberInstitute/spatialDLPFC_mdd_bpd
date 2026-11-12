# Stage 1: Freeze environment inputs and R lockfile

This stage records the package inputs used to construct the MBv container
images. Python's frozen environment files are in [`python/`](python/). The R
lockfile, `renv.lock`, and renv project activation files remain at the project
root:

```text
spatialDLPFC_mdd_bpd/
├── .Rprofile
├── renv.lock
├── renv/
└── container/
    └── 01_environment/
        ├── README.md       # YOU ARE HERE; R update instructions
        └── python/
            └── README.md   # Python update instructions
```

R resources must remain at the project root to function properly. 
The project-level `.Rprofile` activates the project's renv setup, and the
root-level `renv.lock` records the package versions. 

The R Dockerfile and Apptainer definitions copy that lockfile from the project 
root during image builds. Do not move or edit a second copy under this directory;
update the project-root lockfile and rebuild the R image from the project root 
(instructions below).


## Add an R package and update the lockfile

Work from the MBv project root so the project's `.Rprofile` activates renv.
On JHPCE, first load the established R module, then start R:

```bash
cd /path/to/spatialDLPFC_mdd_bpd
module load conda_R/4.4.x
R
```

On another system, start a compatible R session from the project root. In R,
check the current project state, install the package, snapshot the updated
project library, and check the state again:

```r
renv::status()
renv::install("PACKAGE_NAME")
renv::snapshot(prompt = FALSE)
renv::status()
```

If the first `renv::status()` reports an inconsistent project, resolve that
state before installing or snapshotting so unrelated package differences are
not captured as part of the new change.

Replace `PACKAGE_NAME` with the package the code now uses. For example, to add
GEOquery from Bioconductor, use `renv::install("bioc::GEOquery")`. The default
implicit snapshot records packages renv detects as used in project files. Add
the package to the relevant project code before snapshotting; if it is only
installed but not referenced, the default snapshot may not include it. The
`renv::install()` command installs currently available package versions, so
inspect the resulting lockfile diff for dependency upgrades as well as the new
package. See the [renv install reference](https://rstudio.github.io/renv/reference/install.html)
and [snapshot reference](https://rstudio.github.io/renv/reference/snapshot.html).

Quit R, then review the change from the shell:

```bash
git status --short
git diff -- renv.lock
```

Confirm that the new package and its required dependencies are recorded and
that any version changes are intentional. If the package needs system libraries
to compile or run, review the R system-package lists in the Dockerfile and
Apptainer build helper under `container/02_image_build_test/` as well. Commit the
code and lockfile changes together; then rebuild and smoke-test the R image and
rerun the R code check before distributing it.

The lockfile update changes the R image build inputs only. It does not require
rebuilding the Python image unless the Python environment inputs also change.
