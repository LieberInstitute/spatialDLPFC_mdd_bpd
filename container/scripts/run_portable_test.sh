#!/usr/bin/env bash
set -euo pipefail

profile=${1:?Usage: run_portable_test.sh r|python}
case "$profile" in
  r|python) ;;
  *) echo "Profile must be r or python" >&2; exit 2 ;;
esac

command -v apptainer >/dev/null || {
  echo "Apptainer is required. Run this on a clean Linux host with Apptainer installed." >&2
  exit 127
}

repo_root=$(git rev-parse --show-toplevel)
run_id=$(date -u +%Y%m%dT%H%M%SZ)
out="$repo_root/reproducibility_outputs/apptainer/$run_id/${profile}-portable"
mkdir -p "$out"
image="$out/mbv-${profile}.sif"
definition="$repo_root/container/apptainer/${profile}-portable.def"
build_mode=${APPTAINER_BUILD_MODE:-fakeroot}
case "$build_mode" in
  fakeroot) build_args=(--fakeroot) ;;
  root) build_args=() ;;
  *) echo "APPTAINER_BUILD_MODE must be fakeroot or root" >&2; exit 2 ;;
esac

{
  printf 'profile=%s\n' "$profile"
  printf 'run_id_utc=%s\n' "$run_id"
  printf 'repository=%s\n' "$repo_root"
  printf 'git_branch=%s\n' "$(git -C "$repo_root" branch --show-current)"
  printf 'git_commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  printf 'host=%s\n' "$(uname -a)"
  printf 'build_mode=%s\n' "$build_mode"
  apptainer --version
  printf 'build_command=apptainer build'
  printf ' %q' "${build_args[@]}" "$image" "$definition"
  printf '\nexplicit_test_command=apptainer test %q\n' "$image"
} > "$out/run-info.txt"

build_start=$(date +%s)
if ! apptainer build "${build_args[@]}" "$image" "$definition" > "$out/build.log" 2>&1; then
  build_end=$(date +%s)
  printf 'build_elapsed_seconds=%s\nbuild_status=failed\n' "$((build_end-build_start))" >> "$out/run-info.txt"
  cat "$out/build.log" >&2
  exit 1
fi
build_end=$(date +%s)
printf 'build_elapsed_seconds=%s\nbuild_status=passed\n' "$((build_end-build_start))" >> "$out/run-info.txt"

test_start=$(date +%s)
if ! apptainer test "$image" > "$out/test.log" 2>&1; then
  test_end=$(date +%s)
  printf 'explicit_test_elapsed_seconds=%s\nexplicit_test_status=failed\n' "$((test_end-test_start))" >> "$out/run-info.txt"
  cat "$out/test.log" >&2
  exit 1
fi
test_end=$(date +%s)
printf 'explicit_test_elapsed_seconds=%s\nexplicit_test_status=passed\n' "$((test_end-test_start))" >> "$out/run-info.txt"
sha256sum "$image" > "$out/image.sha256"
printf 'image_size_bytes=%s\n' "$(wc -c < "$image" | tr -d ' ')" >> "$out/run-info.txt"
printf 'image_sha256=%s\n' "$(awk '{print $1}' "$out/image.sha256")" >> "$out/run-info.txt"
printf 'result=passed\n' >> "$out/run-info.txt"
printf 'Portable %s build and post-restore test passed. Results: %s\n' "$profile" "$out"
