#!/usr/bin/env python3
"""Run each configured benchmark on every mesh listed in DIRECTORY/all_mesh.txt.

Usage: python3 benchmarking.py DIRECTORY [--output RUN_DIRECTORY] [--timeout SECONDS]
Edit EXECUTABLES, RESULT_HEADER_EXECUTABLE, and PARAMETERS below before running.
Mesh paths in all_mesh.txt can be absolute or relative to DIRECTORY, one per line.
Blank lines and lines starting with # are ignored.

A fresh run directory contains results.csv, runs.csv, parameters.in, and a copy of
all_mesh.txt. All cases run in meshes/ and save generated meshes there, using a
copy of parameters.in; output filenames must be unique to avoid overwriting.
Per-case logs are kept in logs/. Failed cases are recorded and execution continues.
"""

import argparse
import csv
from datetime import datetime
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


# These paths are absolute because they are based on this script's resolved path.
BUILD_DIRECTORY = Path(__file__).resolve().parent / "build"
EXECUTABLES = [
    BUILD_DIRECTORY / "poly_remesh",
    BUILD_DIRECTORY / "poly_remesh_feature_preserving",
    BUILD_DIRECTORY / "sdf_meshing",
    BUILD_DIRECTORY / "offset_meshing",
    # BUILD_DIRECTORY / "octree_hex_mesh",
]
RESULT_HEADER_EXECUTABLE = BUILD_DIRECTORY / "result_header"

# Written to parameters.in in every executable's working directory.
PARAMETERS = {
    "boundary_weight": 1.0,
    "exact_checks": True,
    "grid_size": 100,
    "bbox_target_sizing": 0.02,
    "octree_cube_cells_per_dim": 10,
    "octree_nb_levels": 2,
    "offset": 0.05,
    "verbose": False,
}


def load_meshes(directory):
    meshes = []
    listing = directory / "all_mesh.txt"
    for number, line in enumerate(listing.read_text(encoding="utf-8-sig").splitlines(), 1):
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        mesh = Path(name)
        if not mesh.is_absolute():
            mesh = directory / mesh
        mesh = mesh.resolve()
        if not mesh.is_file():
            raise ValueError(f"{listing}:{number}: mesh does not exist: {mesh}")
        meshes.append(mesh)
    if not meshes:
        raise ValueError(f"No meshes listed in {listing}")
    return meshes


def check_executable(path):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError(f"Executable path must be absolute: {path}")
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"Executable not found or not executable: {path}")
    return path


def write_parameters(filename):
    with filename.open("w", encoding="utf-8") as output:
        for name, value in PARAMETERS.items():
            if isinstance(value, bool):
                value = "true" if value else "false"
            output.write(f"{name} = {value}\n")


def read_header(executable, run_directory):
    process = subprocess.run(
        [str(executable)], cwd=run_directory, capture_output=True,
        text=True, errors="replace", timeout=30, check=False,
    )
    (run_directory / "result_header.stdout.log").write_text(process.stdout, encoding="utf-8")
    (run_directory / "result_header.stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode:
        raise ValueError(f"result_header exited with code {process.returncode}; see its logs")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError("result_header must print exactly one CSV header line")
    header = next(csv.reader(lines, strict=True))
    # The C++ CSV writer may terminate the record with a comma.
    if header and header[-1] == "":
        header.pop()
    if header[:2] != ["Run Type", "Case Name"] or any(not field for field in header):
        raise ValueError("Unexpected result_header output; rebuild the result_header executable")
    return header


def read_result(stdout_file, header, mesh):
    """Find the last matching CSV row, allowing progress messages around it."""
    result = None
    with stdout_file.open(encoding="utf-8", errors="replace") as output:
        for line in output:
            try:
                row = next(csv.reader([line], strict=True))
            except (csv.Error, StopIteration):
                continue
            if row and row[-1] == "":
                row.pop()
            if len(row) == len(header) and row[0] and row[1] == mesh.stem:
                result = row
    if result is None:
        raise ValueError(
            f"No result CSV row for {mesh.stem!r} with {len(header)} columns; "
            "check stdout.log and rebuild result_header if the schema changed"
        )
    return result


def run_case(executable, mesh, directory, header, timeout, stdout_file, stderr_file):
    with stdout_file.open("w") as stdout, stderr_file.open("w") as stderr:
        try:
            process = subprocess.run(
                [str(executable), str(mesh)], cwd=directory,
                stdout=stdout, stderr=stderr, timeout=timeout, check=False,
            )
        except subprocess.TimeoutExpired:
            return None, "", f"Timed out after {timeout:g} seconds"
        except OSError as error:
            return None, "", str(error)
    if process.returncode:
        return None, process.returncode, f"Exited with code {process.returncode}"
    try:
        return read_result(stdout_file, header, mesh), 0, ""
    except ValueError as error:
        return None, 0, str(error)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path, help="Mesh directory containing all_mesh.txt")
    parser.add_argument("--output", type=Path, help="New run directory (must not already exist)")
    parser.add_argument("--timeout", type=float, help="Time limit per case in seconds; default: no limit")
    args = parser.parse_args(argv)
    if args.timeout is not None and not (0 < args.timeout < float("inf")):
        parser.error("--timeout must be finite and positive")

    try:
        directory = args.directory.expanduser().resolve()
        meshes = load_meshes(directory)
        executables = [check_executable(path) for path in EXECUTABLES]
        if not executables:
            raise ValueError("EXECUTABLES is empty")
        header_executable = check_executable(RESULT_HEADER_EXECUTABLE)
        run_directory = (
            args.output.expanduser().resolve() if args.output else
            directory / "benchmark_results" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        )
        run_directory.mkdir(parents=True, exist_ok=False)
        write_parameters(run_directory / "parameters.in")
        shutil.copyfile(directory / "all_mesh.txt", run_directory / "all_mesh.txt")
        header = read_header(header_executable, run_directory)
    except (OSError, ValueError, csv.Error, subprocess.TimeoutExpired) as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        return 2

    print(f"Run directory: {run_directory}", flush=True)
    successes = failures = 0
    total = len(executables) * len(meshes)
    try:
        logs_directory = run_directory / "logs"
        logs_directory.mkdir()
        meshes_directory = run_directory / "meshes"
        meshes_directory.mkdir()
        shutil.copyfile(run_directory / "parameters.in", meshes_directory / "parameters.in")
        with (run_directory / "results.csv").open("w", newline="", encoding="utf-8") as results_file, \
             (run_directory / "runs.csv").open("w", newline="", encoding="utf-8") as runs_file:
            results = csv.writer(results_file)
            runs = csv.writer(runs_file)
            results.writerow(header)
            runs.writerow(["executable", "mesh", "status", "return_code", "seconds", "result_row", "directory", "stdout_log", "stderr_log", "error"])
            results_file.flush()
            runs_file.flush()
            for executable_number, executable in enumerate(executables, 1):
                for mesh_number, mesh in enumerate(meshes, 1):
                    log_prefix = f"{executable_number:02d}_{executable.name}_{mesh_number:05d}"
                    stdout_file = logs_directory / f"{log_prefix}.stdout.log"
                    stderr_file = logs_directory / f"{log_prefix}.stderr.log"
                    print(f"[{successes + failures + 1}/{total}] {executable.name}: {mesh}", flush=True)
                    start = time.monotonic()
                    row, return_code, error = run_case(executable, mesh, meshes_directory, header, args.timeout, stdout_file, stderr_file)
                    elapsed = time.monotonic() - start
                    result_row = ""
                    if row is not None:
                        results.writerow(row)
                        successes += 1
                        result_row = successes + 1  # One-based CSV record number, including the header.
                    else:
                        failures += 1
                        print(f"  FAILED: {error} (log: {stdout_file})", file=sys.stderr, flush=True)
                    runs.writerow([
                        str(executable), str(mesh), "ok" if row is not None else "failed",
                        return_code, f"{elapsed:.3f}", result_row,
                        "meshes", str(stdout_file.relative_to(run_directory)),
                        str(stderr_file.relative_to(run_directory)), error,
                    ])
                    results_file.flush()
                    runs_file.flush()
    except KeyboardInterrupt:
        print(f"\nInterrupted. Completed results are saved in {run_directory}", file=sys.stderr)
        return 130
    except OSError as error:
        print(f"Run failed: {error}. Completed results are saved in {run_directory}", file=sys.stderr)
        return 1

    print(f"Finished: {successes} succeeded, {failures} failed. Results: {run_directory / 'results.csv'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
