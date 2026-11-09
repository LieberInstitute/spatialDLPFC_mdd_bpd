# MBv pipeline container profiles

The project uses separate R and Python environments, so the portable build
uses separate Apptainer images. This also keeps the R image from selecting the
pySCENIC Conda interpreter through `PATH`. The two portable definitions are
`apptainer/r-portable.def` and `apptainer/python-portable.def`; they use the
checked-in `renv.lock`, Conda explicit package lock, and hash-checked pip lock.

These first portable profiles are for Linux x86_64. The Conda explicit lock is
`linux-64`, and the Miniforge bootstrap is the Linux x86_64 installer. The
Docker base tags are currently version tags, not immutable image digests; pin
and record their digests before describing the images as bit-for-bit builds.

## Run on a clean Linux host

Install Apptainer on a disposable/clean Linux host, check out this branch, and
run the profiles separately:

```bash
git switch revision
bash container/scripts/run_portable_test.sh r
bash container/scripts/run_portable_test.sh python
```

The runner defaults to `apptainer build --fakeroot`, which requires fakeroot
support configured in Apptainer. If the clean Linux host has root access and
does not support fakeroot, run the same command with
`APPTAINER_BUILD_MODE=root` from a root shell.

Each run builds an image, runs the definition's `%test` section during the
build, runs `apptainer test` once more, and saves logs and timing metadata under
`reproducibility_outputs/apptainer/<UTC-run-id>/<profile>/`. The output path is
ignored by Git so large images and machine-specific logs stay out of repository
history. Preserve that output directory when reporting build time or a failure.

The build timer includes the build and its `%test` section. The explicit test
time is recorded separately. `run-info.txt` also records the branch, commit,
host, Apptainer version, image size, and SHA-256. Failed build/test logs remain
in place; a failed run exits nonzero.

## Scope and next reproducibility improvements

The tests check the complete R lockfile's package versions and representative
R imports, plus Python, pySCENIC, Arboreto versions/imports and `pip check`.
They do not run the full MBv analysis or validate external study data.

The portable profiles deliberately omit JHPCE-specific APT sandbox and tar
ownership workarounds. If a later JHPCE build needs site-specific behavior, add
it only in a separate JHPCE profile after the portable profiles have been
tested. Do not run `container/scripts/run_jhpce_apptainer_test.sh` as part of
this portable phase.
