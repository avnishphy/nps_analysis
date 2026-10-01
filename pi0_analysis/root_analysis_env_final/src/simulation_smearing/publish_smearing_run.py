#!/usr/bin/env python3
"""Publish validated smearing artifacts with rollback and a committed manifest.

Legacy filenames are individually atomic, not a multi-file read transaction.
Readers needing a consistent set should read the latest manifest once and use
its immutable archive paths. Interrupted legacy publication is recovered before
the next publication (or with --recover-only), under the same advisory locks.
"""

import argparse
from contextlib import ExitStack, contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import sys
import tempfile


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def destination_path(path):
    path = Path(path)
    # Resolve parent aliases without following an existing final-file symlink.
    return path.parent.resolve() / path.name


def validate_destinations(destinations, protected=(), archive=None):
    destinations = [destination_path(path) for path in destinations]
    if len(set(destinations)) != len(destinations):
        raise ValueError("Publication destinations must be distinct, including the manifest")
    protected = {Path(path).resolve() for path in protected}
    for destination in destinations:
        if destination.resolve() in protected:
            raise ValueError(f"Output would replace a protected input/source: {destination}")
        if archive and (archive == destination or archive in destination.parents):
            raise ValueError("Legacy/manifest destinations cannot modify the immutable archive")


def dependency_records(depfiles, extra=()):
    paths = set(Path(path).resolve() for path in extra)
    for depfile in depfiles:
        words = shlex.split(Path(depfile).read_text().replace("\\\n", " "))
        paths.update(Path(word).resolve() for word in words[1:])
    return [{"path": str(path), "sha256": sha256(path)} for path in sorted(paths)]


def dependency_identity(depfiles, extra=()):
    payload = json.dumps(dependency_records(depfiles, extra), sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        if os.path.lexists(temp):
            os.unlink(temp)


def replace_artifact(source, destination):
    """Copy to a unique same-directory temporary, then atomically replace."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix="." + destination.name + ".", dir=destination.parent)
    os.close(fd)
    try:
        if Path(source).is_symlink():
            os.unlink(temp)
            os.symlink(os.readlink(source), temp)
        else:
            shutil.copy2(source, temp)
            with open(temp, "rb") as stream:
                os.fsync(stream.fileno())
        os.replace(temp, destination)
    finally:
        if os.path.lexists(temp):
            os.unlink(temp)


@contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield stream


@contextmanager
def destination_locks(destinations, owner):
    # Sorting gives overlapping custom output sets a consistent lock order.
    with ExitStack() as stack:
        streams = []
        for destination in sorted(set(map(str, destinations))):
            path = Path(destination)
            stream = stack.enter_context(locked(path.parent / ("." + path.name + ".smearing.lock")))
            stream.seek(0)
            previous_owner = stream.read().strip()
            if previous_owner and previous_owner != owner:
                raise RuntimeError(f"Destination {path} belongs to manifest {previous_owner}; "
                                   "refusing cross-manifest publication/recovery")
            streams.append(stream)
        for stream in streams:
            stream.seek(0)
            stream.truncate()
            stream.write(owner + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        yield


def rollback(journal):
    failures = []
    for item in reversed(journal["files"]):
        try:
            if item["existed"]:
                replace_artifact(item["backup"], item["destination"])
            elif os.path.lexists(item["destination"]):
                # This slot was absent before this transaction; remove our file.
                os.unlink(item["destination"])
        except OSError as error:
            failures.append(f"{item['destination']}: {error}")
    if failures:
        raise RuntimeError("Rollback incomplete; retain journal/backups: " + "; ".join(failures))


def recover_pending(recovery_root):
    if not recovery_root.exists():
        return
    for directory in sorted(recovery_root.iterdir()):
        journal_path = directory / "journal.json"
        if not journal_path.exists():
            # Backups may be incomplete, but no destination was touched yet.
            shutil.rmtree(directory)
            continue
        journal = json.loads(journal_path.read_text())
        with destination_locks((item["destination"] for item in journal["files"]),
                               journal["manifest"]):
            if journal["state"] != "committed":
                rollback(journal)
                print(f"[recovered] Previous legacy outputs from {directory}")
            shutil.rmtree(directory)


def publish(archive_staging, archive, manifest, pairs, run_id, protected=()):
    archive_staging, archive, manifest = map(Path, (archive_staging, archive, manifest))
    archive_staging = archive_staging.resolve()
    archive = destination_path(archive)
    manifest = destination_path(manifest)
    # Sources must be artifacts already assembled into this validated archive.
    artifacts = []
    for source, destination in pairs:
        source = Path(source).resolve()
        relative = source.relative_to(archive_staging)
        if not source.is_file():
            raise ValueError(f"Missing publication artifact: {source}")
        artifacts.append({"archive_path": str(archive / relative),
                          "relative_path": str(relative),
                          "destination": str(destination_path(destination)),
                          "sha256": sha256(source)})
    destinations = [item["destination"] for item in artifacts] + [str(manifest)]
    validate_destinations(destinations, protected, archive)
    recovery_root = manifest.parent / ".smearing_publication_recovery"
    with locked(manifest.parent / ".smearing_publication.lock"):
        recover_pending(recovery_root)
        if archive.exists():
            raise FileExistsError(f"Refusing to replace immutable run archive: {archive}")
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive_staging.rename(archive)
        print(f"[archive committed] {archive}", flush=True)
        # Everything from here may fail without removing the complete archive.
        with destination_locks(destinations, str(manifest)):
            recovery_root.mkdir(parents=True, exist_ok=True)
            recovery = Path(tempfile.mkdtemp(prefix="publication.", dir=recovery_root))
            journal = {"state": "prepared", "archive": str(archive),
                       "manifest": str(manifest), "files": []}
            try:
                for index, destination in enumerate(destinations):
                    path = Path(destination)
                    existed = os.path.lexists(path)
                    if existed and not path.is_file() and not path.is_symlink():
                        raise ValueError(f"Output destination is not a regular file/symlink: {path}")
                    backup = recovery / str(index)
                    if existed:
                        shutil.copy2(path, backup, follow_symlinks=False)
                    journal["files"].append({"destination": destination,
                                             "existed": existed, "backup": str(backup)})
                atomic_json(recovery / "journal.json", journal)
            except BaseException:
                shutil.rmtree(recovery)
                raise
            try:
                for item in artifacts:
                    replace_artifact(item["archive_path"], item["destination"])
                    if sha256(item["destination"]) != item["sha256"]:
                        raise RuntimeError(f"Published checksum mismatch: {item['destination']}")
                atomic_json(manifest, {"version": 1, "run_id": run_id,
                                       "archive": str(archive), "artifacts": artifacts})
                journal["state"] = "committed"
                atomic_json(recovery / "journal.json", journal)
            except BaseException:
                # SIGKILL cannot be caught; its on-disk journal remains recoverable.
                rollback(journal)
                shutil.rmtree(recovery)
                raise
            shutil.rmtree(recovery)
        print(f"[latest committed] {manifest}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-staging")
    parser.add_argument("--archive")
    parser.add_argument("--manifest")
    parser.add_argument("--run-id")
    parser.add_argument("--file", nargs=2, action="append", default=[])
    parser.add_argument("--recover-only", action="store_true")
    parser.add_argument("--check-destinations", action="store_true")
    parser.add_argument("--destination", action="append", default=[])
    parser.add_argument("--protected-path", action="append", default=[])
    parser.add_argument("--dependency-identity", action="store_true")
    parser.add_argument("--depfile", action="append", default=[])
    parser.add_argument("--extra-source", action="append", default=[])
    parser.add_argument("--build-provenance")
    parser.add_argument("--executable", action="append", default=[])
    parser.add_argument("--toolchain-info")
    args = parser.parse_args()
    if args.check_destinations:
        dependencies = [item["path"] for item in dependency_records(args.depfile, args.extra_source)]
        validate_destinations(args.destination, args.protected_path + dependencies)
        return 0
    if args.dependency_identity:
        print(dependency_identity(args.depfile, args.extra_source))
        return 0
    if args.build_provenance:
        atomic_json(args.build_provenance, {
            "dependencies": dependency_records(args.depfile, args.extra_source),
            "executables": [{"path": str(Path(path).resolve()), "sha256": sha256(path)}
                            for path in args.executable],
            "toolchain": Path(args.toolchain_info).read_text()})
        return 0
    if not args.manifest:
        parser.error("--manifest is required")
    if args.recover_only:
        manifest = destination_path(args.manifest)
        with locked(manifest.parent / ".smearing_publication.lock"):
            recover_pending(manifest.parent / ".smearing_publication_recovery")
        return 0
    if not all((args.archive_staging, args.archive, args.run_id, args.file)):
        parser.error("publication requires archive-staging, archive, run-id and file entries")
    dependencies = [item["path"] for item in dependency_records(args.depfile, args.extra_source)]
    publish(args.archive_staging, args.archive, args.manifest, args.file, args.run_id,
            args.protected_path + dependencies)
    return 0


if __name__ == "__main__":
    def interrupted(signum, _frame):
        raise KeyboardInterrupt(f"Interrupted by signal {signum}")
    signal.signal(signal.SIGTERM, interrupted)
    try:
        sys.exit(main())
    except (Exception, KeyboardInterrupt) as error:
        print(f"[ERROR] Publication/provenance failed: {error}", file=sys.stderr)
        sys.exit(1)
