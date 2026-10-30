#!/usr/bin/env python3
"""Create a hash-checked pip lock from the MBv export's exact version pins.

The source export already records every installed pip package, so this script
does not resolve dependency versions. It records the SHA-256 digests published
for each exact PyPI release and adds the local pySCENIC/Arboreto wheels.
"""

from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "requirements-pip.in"
OUTPUT = ROOT / "requirements-pip.lock"


def parse_pin(line: str) -> tuple[str, str]:
    name, separator, version = line.strip().partition("==")
    if not separator or not name or not version:
        raise ValueError(f"Expected an exact name==version pin, got: {line!r}")
    return name, version


def pypi_project_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def release_hashes(pin: tuple[str, str]) -> tuple[str, str, list[str]]:
    name, version = pin
    url = f"https://pypi.org/pypi/{pypi_project_name(name)}/{version}/json"
    request = Request(url, headers={"User-Agent": "MBv-pipeline-lock/1.0"})
    with urlopen(request, timeout=45) as response:
        release = json.load(response)

    wheels = sorted(
        {item["digests"]["sha256"] for item in release["urls"]
         if item["filename"].endswith(".whl") and item.get("digests", {}).get("sha256")}
    )
    artifacts = wheels or sorted(
        {item["digests"]["sha256"] for item in release["urls"]
         if item.get("digests", {}).get("sha256")}
    )
    if not artifacts:
        raise RuntimeError(f"PyPI returned no files with SHA-256 digests for {name}=={version}")
    return name, version, artifacts


def requirement_block(requirement: str, hashes: list[str]) -> str:
    first = f"{requirement} \\" if hashes else requirement
    if not hashes:
        return first
    return first + "\n" + "\n".join(f"    --hash=sha256:{digest}" for digest in hashes)


def main() -> None:
    pins = [line.strip() for line in INPUT.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    parsed = [parse_pin(line) for line in pins]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(release_hashes, parsed))

    blocks = ["# Generated from requirements-pip.in using PyPI release SHA-256 digests.",
              "# Install from this directory with: python -m pip install --no-deps -r requirements-pip.lock",
              "--require-hashes"]
    for name, version, hashes in results:
        blocks.append(requirement_block(f"{name}=={version}", hashes))

    for wheel in sorted((ROOT / "wheels").glob("*.whl")):
        digest = sha256(wheel.read_bytes()).hexdigest()
        blocks.append(requirement_block(f"./wheels/{wheel.name}", [digest]))

    OUTPUT.write_text("\n\n".join(blocks) + "\n")
    print(f"Wrote {OUTPUT} ({len(parsed)} PyPI pins, {len(list((ROOT / 'wheels').glob('*.whl')))} local wheels)")


if __name__ == "__main__":
    main()
