from importlib.metadata import version
import platform

expected = {
    "python": "3.10.14",
    "pyscenic": "0.12.1+11.g06bafba",
    "arboreto": "0.1.6",
}
actual = {
    "python": platform.python_version(),
    "pyscenic": version("pyscenic"),
    "arboreto": version("arboreto"),
}
for package, expected_version in expected.items():
    print(f"{package}: expected={expected_version} actual={actual[package]}")
    if actual[package] != expected_version:
        raise SystemExit(f"Version mismatch for {package}")

import arboreto.algo  # noqa: F401
import dask  # noqa: F401
import pyscenic  # noqa: F401
print("Python environment imports passed")
