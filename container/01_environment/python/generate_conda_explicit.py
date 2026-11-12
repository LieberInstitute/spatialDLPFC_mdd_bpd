#!/usr/bin/env python3
"""Render the exact Conda builds from the exported environment spec.

conda-lock can add transitive Conda packages that were not explicit in the
export. This script keeps only the Conda packages listed in
environment-conda.yml and borrows their artifact URLs and MD5 checksums from
conda-lock-conda.yml. Pip packages are installed separately from their hash lock.
"""

from pathlib import Path
from urllib.parse import urlparse

import yaml


ROOT = Path(__file__).resolve().parent
ENVIRONMENT = ROOT / "environment-conda.yml"
LOCK = ROOT / "conda-lock-conda.yml"
OUTPUT = ROOT / "conda-linux-64.lock"


def main() -> None:
    environment = yaml.safe_load(ENVIRONMENT.read_text())
    lock = yaml.safe_load(LOCK.read_text())
    locked_packages = {
        package["name"]: package
        for package in lock["package"]
        if package["manager"] == "conda" and package["platform"] == "linux-64"
    }

    lines = [
        "# Linux x86_64 explicit lock for the exact Conda builds in environment-conda.yml.",
        "# Package URLs and MD5 checksums are taken from conda-lock-conda.yml.",
        "# Create with: conda create --name pyscenic_bioconda --file conda-linux-64.lock",
        "@EXPLICIT",
    ]
    count = 0
    for dependency in environment["dependencies"]:
        if not isinstance(dependency, str):
            continue
        parts = dependency.split("=")
        if len(parts) != 3 or not all(parts):
            raise ValueError(f"Expected an exact Conda version/build pin, got {dependency!r}")
        name, version, build = parts
        package = locked_packages.get(name)
        if package is None or package["version"] != version:
            raise ValueError(f"No matching locked Conda package for {dependency!r}")
        url = package["url"]
        filename = Path(urlparse(url).path).name
        if not filename.startswith(f"{name}-{version}-{build}."):
            raise ValueError(f"Locked artifact does not match {dependency!r}: {filename}")
        md5 = package["hash"].get("md5")
        if not md5:
            raise ValueError(f"No MD5 checksum recorded for {dependency!r}")
        lines.append(f"{url}#{md5}")
        count += 1

    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUTPUT} with {count} exact Conda package artifacts")


if __name__ == "__main__":
    main()
