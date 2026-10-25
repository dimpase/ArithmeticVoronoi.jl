"""Restore frozen D5/D8 regression fixtures from an explicit data checkout.

Existing files with different contents are never overwritten. Use --check
to verify an already restored directory without writing anything.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

PROJECT = Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def restore(data_root, directory, check_only=False):
    manifest_path = data_root / "archives/fixtures/two_cells_v1.json"
    manifest = json.loads(manifest_path.read_text())
    archive = manifest_path.parent / manifest["archive"]
    if not archive.is_file() or digest(archive.read_bytes()) != manifest["sha256"]:
        raise ValueError("Fixture archive is missing or has the wrong hash. Run git lfs pull in the selected ArithmeticVoronoiData.jl checkout.")
    directory = directory.resolve()
    with tarfile.open(archive) as source:
        members = source.getmembers()
        if len(members) != len(manifest["files"]) or {m.name for m in members} != set(manifest["files"]):
            raise ValueError("Fixture archive members do not match the manifest")
        missing = []
        for member in members:
            destination = (directory / member.name).resolve()
            if not member.isfile() or not destination.is_relative_to(directory):
                raise ValueError(f"Invalid fixture archive member: {member.name}")
            expected = manifest["files"][member.name]
            if digest(source.extractfile(member).read()) != expected:
                raise ValueError(f"Fixture checksum mismatch: {member.name}")
            if destination.exists():
                if not destination.is_file() or digest(destination.read_bytes()) != expected:
                    raise ValueError(f"Existing fixture differs; refusing to overwrite: {destination}")
            else:
                missing.append(member)
        if missing and check_only:
            raise ValueError(f"Missing {len(missing)} fixtures; run scripts/restore_fixtures.py")
        if missing:
            source.extractall(directory, members=missing, filter="data")
    print(f"Verified {len(members)} fixtures; restored {len(missing)} files to {directory}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True,
                        help="checkout of https://github.com/dimpase/ArithmeticVoronoiData.jl")
    parser.add_argument("--directory", type=Path, default=PROJECT / "results/two_cells")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        restore(args.data_root, args.directory, args.check)
    except (OSError, ValueError, tarfile.TarError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
