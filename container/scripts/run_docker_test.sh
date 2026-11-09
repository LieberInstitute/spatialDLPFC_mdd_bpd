#!/usr/bin/env bash
set -euo pipefail

profile=${1:?Usage: run_docker_test.sh r|python}
case "$profile" in
  r|python) ;;
  *) echo "Profile must be r or python" >&2; exit 2 ;;
esac

if [[ -n "${DOCKER_BIN:-}" ]]; then
  docker_cmd=("$DOCKER_BIN")
elif command -v docker >/dev/null 2>&1; then
  docker_cmd=(docker)
elif [[ -x /Applications/Docker.app/Contents/Resources/bin/docker ]]; then
  docker_cmd=(/Applications/Docker.app/Contents/Resources/bin/docker)
  export PATH="/Applications/Docker.app/Contents/Resources/bin:$PATH"
else
  echo "Docker CLI not found. Install/start Docker Desktop on macOS or Docker Engine on Linux." >&2
  exit 127
fi

"${docker_cmd[@]}" info >/dev/null || {
  echo "Docker daemon is unavailable. Start Docker Desktop/Engine and retry." >&2
  exit 1
}

repo_root=$(git rev-parse --show-toplevel)
run_id=$(date -u +%Y%m%dT%H%M%SZ)
out="$repo_root/reproducibility_outputs/docker/$run_id/$profile"
mkdir -p "$out"
platform=linux/amd64
tag="mbv-${profile}:test"
dockerfile="$repo_root/container/docker/$profile/Dockerfile"
build_cmd=("${docker_cmd[@]}" build --pull --platform "$platform" \
  -f "$dockerfile" -t "$tag" "$repo_root")

{
  printf 'profile=%s\nrun_id_utc=%s\n' "$profile" "$run_id"
  printf 'repository=%s\n' "$repo_root"
  printf 'git_branch=%s\n' "$(git -C "$repo_root" branch --show-current)"
  printf 'git_commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  printf 'host=%s\nplatform=%s\n' "$(uname -a)" "$platform"
  "${docker_cmd[@]}" version --format 'Docker client={{.Client.Version}} server={{.Server.Version}}'
  printf 'build_command='
  printf '%q ' "${build_cmd[@]}"
  printf '\n'
} > "$out/run-info.txt"

build_start=$(date +%s)
if ! "${build_cmd[@]}" > "$out/build.log" 2>&1; then
  build_end=$(date +%s)
  printf 'build_elapsed_seconds=%s\nbuild_status=failed\n' "$((build_end-build_start))" >> "$out/run-info.txt"
  cat "$out/build.log" >&2
  exit 1
fi
build_end=$(date +%s)
printf 'build_elapsed_seconds=%s\nbuild_status=passed\n' "$((build_end-build_start))" >> "$out/run-info.txt"

test_start=$(date +%s)
if [[ "$profile" == r ]]; then
  if ! "${docker_cmd[@]}" run --rm --platform "$platform" "$tag" Rscript --vanilla \
      /opt/mbv/scripts/verify_r_lockfile.R /opt/mbv/renv-library /opt/mbv/project/renv.lock \
      > "$out/test.log" 2>&1 || \
      ! "${docker_cmd[@]}" run --rm --platform "$platform" "$tag" Rscript --vanilla -e \
      'library(sf); library(s2); library(Seurat); library(SingleCellExperiment); library(spatialLIBD)' \
      >> "$out/test.log" 2>&1; then
    test_end=$(date +%s)
    printf 'smoke_test_elapsed_seconds=%s\nsmoke_test_status=failed\n' "$((test_end-test_start))" >> "$out/run-info.txt"
    cat "$out/test.log" >&2
    exit 1
  fi
else
  if ! "${docker_cmd[@]}" run --rm --platform "$platform" "$tag" sh -eu -c \
      'python /opt/mbv/scripts/verify_python_env.py && python -m pip check' \
      > "$out/test.log" 2>&1; then
    test_end=$(date +%s)
    printf 'smoke_test_elapsed_seconds=%s\nsmoke_test_status=failed\n' "$((test_end-test_start))" >> "$out/run-info.txt"
    cat "$out/test.log" >&2
    exit 1
  fi
fi
test_end=$(date +%s)
printf 'smoke_test_elapsed_seconds=%s\nsmoke_test_status=passed\n' "$((test_end-test_start))" >> "$out/run-info.txt"

image_id=$("${docker_cmd[@]}" image inspect "$tag" --format '{{.Id}}')
image_size=$("${docker_cmd[@]}" image inspect "$tag" --format '{{.Size}}')
printf 'image_id=%s\nimage_size_bytes=%s\nresult=passed\n' \
  "$image_id" "$image_size" >> "$out/run-info.txt"
printf 'Docker %s build and smoke test passed. Results: %s\n' "$profile" "$out"
