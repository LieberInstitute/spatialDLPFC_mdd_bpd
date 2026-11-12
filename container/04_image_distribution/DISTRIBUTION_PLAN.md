# MBv container distribution plan

## Purpose and scope

The MBv pipeline currently uses two independent software environments for the
bulk of its analysis in code directories `00`–`10`:

- **R image:** R 4.4.3 and the package versions recorded in `renv.lock`.
- **Python image:** the locked `pyscenic_bioconda` environment, including
  pySCENIC and Arboreto, with the pySCENIC Dask memory patch retained.

The tensorQTL work in `11_eQTL_coloc` is outside the current image scope. It may
need a separate image later. See [container build guide](../README.md)
for image build instructions, supported platforms, test status, and measured
build resources.

## Distribution options

| Option | What users receive | Best fit | Considerations |
| --- | --- | --- | --- |
| Build from the project checkout | Docker image or Apptainer SIF built locally from the checked-in definitions and locks | Users who can build containers and want to inspect or modify the build | Requires a project clone, network access during the build, substantial build time, and enough memory and disk. |
| Publish Docker images to GHCR | Versioned OCI/Docker images, pulled with Docker | Workstations using Docker Desktop or Docker Engine, including macOS and Windows | Convenient for users; public images can be pulled without a GitHub login. The Docker images target `linux/amd64`; Apple Silicon runs that target using emulation. |
| Publish Apptainer images as SIF files | Ready-to-run `.sif` files, for example through GitHub Releases or an institutional storage location | Linux and HPC users who need Apptainer | SIF files are built and tested separately for portable Linux and JHPCE. Publishing a SIF through a registry requires an explicit OCI/ORAS distribution and retrieval procedure; it is not the same as publishing a Docker image. |
| Share image files through Globus | `.sif` files directly, and Docker images exported as `.tar.gz` archives | Project collaborators who already access MBv data through Globus, especially HPC users | Globus transfers files from a configured collection; users must download the artifacts before running them. It is not a Docker or Apptainer registry, and collection sharing permissions and subscription eligibility must be confirmed. |
| Publish Apptainer-compatible OCI images | OCI images that Apptainer can pull and convert to SIF | Linux users who prefer registry-based delivery | Must be tested with the selected registry, architecture, tags, and Apptainer version. This is a later option; the existing Docker images should not be assumed to be the tested Apptainer builds. |

## Recommendation for this project

Use **GHCR as the primary distribution route for the two Docker images**. Publish
the R and Python images as separate packages, under stable, lowercase names:

```text
ghcr.io/lieberinstitute/mbv-r
ghcr.io/lieberinstitute/mbv-python
```

Publish only images that passed their corresponding final code and environment
checks. Use a release tag such as `v1.0.0`, and also record the full MBv Git
commit associated with each image. Do not use `latest` as the reproducibility
reference: tags can be moved, while a digest identifies the exact image content.

### GHCR layer-size check

GitHub documents a **10 GB maximum size per image layer** and a 10-minute
upload timeout; the documented size limit is per layer, not a 10 GB cap on the
whole image ([GHCR documentation](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)).
The current Docker image totals are about 6.80 GB for R and 13.60 GB for Python.
Local `docker image history` output reports the largest layer as 2.73 GB for R
and 9.08 GB for Python. On those measurements, neither image exceeds the layer
limit, although Python is closer to it. These sizes are from local Docker
history, not a completed GHCR upload; verify the actual push and watch for
upload timeouts before announcing a release. Recheck layer sizes after any
Dockerfile or dependency changes.

### Globus as an alternative file-transfer route

Because the project already uses Globus to share large data files, it may also
be used to make container artifacts available to collaborators. This is a
file-download workflow: Globus transfers files between storage collections,
but does not provide `docker pull` or `apptainer pull` access to a shared folder.
Users need the appropriate container runtime installed after downloading the
artifacts.

- **Apptainer:** share `.sif` files and a share-ready `SHA256SUMS` file in a
  read-only collection folder. Users download the SIFs and checksum file,
  verify them with `sha256sum -c`, then run the SIFs with Apptainer. The current
  build runner writes the original absolute image path in `image.sha256`; make
  a share-ready checksum file using relative filenames before transfer, for
  example from the UTC run directory:

  ```bash
  sha256sum r-jhpce/mbv-r-jhpce.sif \
    python-jhpce/mbv-python-jhpce.sif > SHA256SUMS
  ```

  The latest JHPCE SIFs are about 1.73 GiB (R) and 3.79 GiB (Python); a
  portable profile should be offered separately after it is tested on clean
  Linux.
- **Docker:** export the versioned `linux/amd64` images to one archive, then
  transfer that archive and its checksum through Globus. For the current local
  test tags, the export command is:

  ```bash
  docker image save --platform linux/amd64 mbv-r:test mbv-python:test \
    | gzip > mbv-docker-images.tar.gz
  sha256sum mbv-docker-images.tar.gz > mbv-docker-images.tar.gz.sha256
  ```

  After downloading, a user checks the checksum and loads the archive with
  `docker load --input mbv-docker-images.tar.gz`. Docker supports saving images
  to tar archives and loading gzip-compressed archives ([save](https://docs.docker.com/reference/cli/docker/image/save/),
  [load](https://docs.docker.com/reference/cli/docker/image/load/)). The
  compressed archive size has not yet been measured; users need enough local
  disk space for the download and loaded images.

Globus guest collection sharing is subject to the source collection's
configuration and subscription eligibility. Globus documents guest sharing as
available on collections with an active subscription, with some Globus Connect
Personal cases requiring subscription-group membership
([Globus sharing guide](https://docs.globus.org/guides/tutorials/manage-files/share-files/)).
Check with the collection administrator that a read-only software-artifact
folder can be shared with the intended user group. Existing data-sharing access
does not by itself confirm that the same paths or collection permissions can
be used for container images.

For this project, Globus is a useful **secondary delivery option** for
collaborators who already use the project collection or need a direct SIF
download. GHCR remains simpler for workstation users because Docker pulls
versioned images directly, and Apptainer can use a registry route when that is
validated. Document artifact versions, source commit, profile, checksum, and
required runtime with any Globus-shared image files.

Keep Docker and Apptainer distribution decisions separate. GHCR can store Docker
and OCI images; it does not make the tested JHPCE SIF available automatically.
For Apptainer, retain the current build-and-test process until a clean Linux
portable build has been validated. After that, choose between publishing SIF
files as release assets/institutional files or publishing a separately validated
OCI image for Apptainer to pull. Document each artifact's profile (`portable` or
`jhpce`) and test result so users do not mistake one for the other.

The Dockerfiles currently identify the image by title, version, and description.
Before a first GHCR release, consider adding the OCI label
`org.opencontainers.image.source=https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd`
to associate each package with the source repository. Because changing a label
changes the built image, rebuild and rerun the relevant checks before publishing
that labeled image. Alternatively, publish the already-tested image and connect
the package to the repository in GitHub package settings.

## Publishing Docker images to GHCR

These instructions publish the **already-built local images** named
`mbv-r:test` and `mbv-python:test`. First confirm that these tags still refer to
the images that passed the final checks, and note the exact project commit:

```bash
git status --short
git rev-parse HEAD
docker image inspect mbv-r:test --format '{{.Id}} {{.Os}}/{{.Architecture}}'
docker image inspect mbv-python:test --format '{{.Id}} {{.Os}}/{{.Architecture}}'
```

The expected platform is `linux/amd64`. If there are uncommitted changes or the
image tags do not point to the tested builds, stop and build/test the intended
images before publishing.

### Authenticate

Create a GitHub personal access token (classic) with the `write:packages` scope.
If the Lieber Institute organization requires it, authorize the token for that
organization's SSO. Sign in through Docker's prompt so the token is not embedded
in shell history:

```bash
docker login ghcr.io -u YOUR_GITHUB_USERNAME
```

Enter the token at the password prompt. The account must also have permission to
publish packages for the organization. GitHub's current instructions are in
[Working with the Container registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry).

### Tag and push

Choose a new release tag for the tested pair of images. Replace the example
`v1.0.0` with the project's actual release version. The optional second tag
records the source commit in the registry; use the short commit hash obtained
above. Never reuse a release tag for different image content.

```bash
VERSION=v1.0.0
COMMIT=$(git rev-parse --short=12 HEAD)

docker tag mbv-r:test \
  "ghcr.io/lieberinstitute/mbv-r:${VERSION}"
docker tag mbv-r:test \
  "ghcr.io/lieberinstitute/mbv-r:${COMMIT}"

docker tag mbv-python:test \
  "ghcr.io/lieberinstitute/mbv-python:${VERSION}"
docker tag mbv-python:test \
  "ghcr.io/lieberinstitute/mbv-python:${COMMIT}"

docker push "ghcr.io/lieberinstitute/mbv-r:${VERSION}"
docker push "ghcr.io/lieberinstitute/mbv-r:${COMMIT}"
docker push "ghcr.io/lieberinstitute/mbv-python:${VERSION}"
docker push "ghcr.io/lieberinstitute/mbv-python:${COMMIT}"
```

The tags make the local image available under the GHCR name; `docker push`
uploads it. Publishing creates/updates packages under the organization, so
confirm the names and release tag before pushing.

### Set visibility and verify

After the first push, inspect each package under the Lieber Institute's GitHub
**Packages** page. Set it to **Public** if the intended users should be able to
pull anonymously. GHCR public container packages can be pulled without signing
in. See [GitHub package visibility and access](https://docs.github.com/en/packages/learn-github-packages/configuring-a-packages-access-control-and-visibility).

Pull the published release tags to verify registry access:

```bash
docker pull "ghcr.io/lieberinstitute/mbv-r:${VERSION}"
docker pull "ghcr.io/lieberinstitute/mbv-python:${VERSION}"
```

Record the release tag, source commit, package URLs, and image digests with the
release notes. Get the digest from the registry with:

```bash
docker buildx imagetools inspect "ghcr.io/lieberinstitute/mbv-r:${VERSION}"
docker buildx imagetools inspect "ghcr.io/lieberinstitute/mbv-python:${VERSION}"
```

Users should use versioned image references (or the recorded digest), not a
floating tag, when they need to reproduce a run.

## Pulling and using a published Docker image

After the packages are public, a user can pull them without a GitHub login:

```bash
docker pull ghcr.io/lieberinstitute/mbv-r:v1.0.0
docker pull ghcr.io/lieberinstitute/mbv-python:v1.0.0
```

To run a project script, clone the project at the matching code revision and
mount the working directory into the container. For example, from the project
root:

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  ghcr.io/lieberinstitute/mbv-r:v1.0.0 \
  Rscript --vanilla path/to/script.R
```

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" \
  -w /work \
  ghcr.io/lieberinstitute/mbv-python:v1.0.0 \
  python3 -u path/to/script.py
```

Replace the example version and script paths with the published release and the
script to run. The image supplies the software environment; the project code
and any input data still need to be available on the host and mounted as needed.
Users must clone the project at the matching commit for the code and image to
correspond. These commands illustrate container invocation only; they do not
define or execute the full MBv workflow.

## Future Nextflow use

Nextflow can refer to the published R and Python image URIs in process-level
container settings. Each process can select the image that matches the language
and environment it needs. Pin the image to a release tag or, preferably for
strict reproducibility, its digest. The workflow implementation and process
configuration are out of scope for this distribution plan; the separate-image
approach already supports that later integration.

## Release checklist

- [ ] Confirm the local R and Python image tags correspond to the final tested
      builds and record the source commit.
- [ ] Choose a new version tag; do not overwrite an existing release tag.
- [ ] Authenticate to GHCR with an account authorized to publish for the
      Lieber Institute organization.
- [ ] Push both images and set their intended visibility.
- [ ] Verify pulls and record image digests alongside the code commit.
- [ ] Update project documentation with package links, released tags, digests,
      and any known platform limitations.
- [ ] Keep Apptainer publication as a separate release decision until portable
      Linux testing and the delivery method have been completed.
