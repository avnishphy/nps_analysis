"""Bounded publication failure/recovery checks; no ROOT or production inputs."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
HELPER = REPO / 'src/simulation_smearing/publish_smearing_run.py'
spec = importlib.util.spec_from_file_location('publish', HELPER)
pub = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pub)


class PublicationChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='publication_test_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.stage = self.root / 'stage'
        self.stage.mkdir()
        self.archive = self.root / 'archives/run1'
        self.output = self.root / 'latest'
        self.output.mkdir()
        self.manifest = self.output / 'smearing_latest.json'
        self.manifest.write_text('{"run_id":"old"}\n')
        self.pairs = []
        for index in range(3):
            source = self.stage / f'{index}.root'
            source.write_bytes(f'new-{index}'.encode())
            destination = self.output / source.name
            if index != 1:
                destination.write_bytes(f'old-{index}'.encode())
            self.pairs.append((source, destination))

    def assert_previous(self):
        self.assertEqual(self.pairs[0][1].read_bytes(), b'old-0')
        self.assertFalse(self.pairs[1][1].exists())
        self.assertEqual(self.pairs[2][1].read_bytes(), b'old-2')
        self.assertEqual(self.manifest.read_text(), '{"run_id":"old"}\n')
        self.assertEqual((self.archive/'0.root').read_bytes(), b'new-0')
        self.assertEqual((self.archive/'1.root').read_bytes(), b'new-1')
        self.assertEqual((self.archive/'2.root').read_bytes(), b'new-2')

    def test_success(self):
        pub.publish(self.stage, self.archive, self.manifest, self.pairs, 'run1')
        manifest = json.loads(self.manifest.read_text())
        self.assertEqual(manifest['run_id'], 'run1')
        for item in manifest['artifacts']:
            self.assertEqual(pub.sha256(item['archive_path']), item['sha256'])
            self.assertEqual(pub.sha256(item['destination']), item['sha256'])
        self.assertFalse(self.stage.exists())

    def test_mid_publish_failure_rolls_back_and_retains_archive(self):
        original = pub.replace_artifact
        calls = 0
        def fail_after_second(source, destination):
            nonlocal calls
            calls += 1
            original(source, destination)
            if calls == 2:
                raise OSError('injected copy/rename failure after second replacement')
        with patch.object(pub, 'replace_artifact', fail_after_second):
            with self.assertRaises(OSError):
                pub.publish(self.stage, self.archive, self.manifest, self.pairs, 'run1')
        self.assert_previous()

    def test_interruption_journal_recovery(self):
        program = '''
import importlib.util, os, sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('publish',sys.argv[1])
p=importlib.util.module_from_spec(spec); spec.loader.exec_module(p)
r=Path(sys.argv[2]); original=p.replace_artifact
def interrupted(s,d):
    original(s,d)
    os._exit(77)
p.replace_artifact=interrupted
p.publish(r/'stage',r/'archives/run1',r/'latest/smearing_latest.json',
          [(r/'stage'/f'{i}.root',r/'latest'/f'{i}.root') for i in range(3)],'run1')
'''
        result = subprocess.run([sys.executable, '-c', program, str(HELPER), str(self.root)])
        self.assertEqual(result.returncode, 77)
        self.assertEqual(self.pairs[0][1].read_bytes(), b'new-0')
        self.assertEqual(self.manifest.read_text(), '{"run_id":"old"}\n')
        subprocess.run([sys.executable, str(HELPER), '--manifest', str(self.manifest),
                        '--recover-only'], check=True)
        self.assert_previous()

    def test_same_size_timestamp_preserving_change_is_detected(self):
        path = self.stage / 'small.conf'
        path.write_bytes(b'cut=1\n')
        first_stat = path.stat()
        script = (REPO/'src/simulation_smearing/run_smearing_pipeline.sh').read_text()
        start = script.index('file_identity() {')
        function = script[start:script.index('\n}\n', start)+3]
        command = 'set -o pipefail\n' + function + '\nfile_identity "$1"'
        def identity():
            return subprocess.check_output(['bash', '-c', command, '--', str(path)])
        first = identity()
        path.write_bytes(b'cut=2\n')
        os.utime(path, ns=(first_stat.st_atime_ns, first_stat.st_mtime_ns))
        self.assertEqual(path.stat().st_size, first_stat.st_size)
        self.assertEqual(path.stat().st_mtime_ns, first_stat.st_mtime_ns)
        self.assertNotEqual(identity(), first)

    def test_shell_identity_missing_input_propagates_failure(self):
        script = (REPO/'src/simulation_smearing/run_smearing_pipeline.sh').read_text()
        start = script.index('file_identity() {')
        function = script[start:script.index('\n}\n', start)+3]
        command = ('set -o pipefail\n' + function +
                   '\nidentity="$(file_identity "$1")"; rc=$?; printf "%s" "$identity"; exit "$rc"')
        result = subprocess.run(['bash', '-c', command, '--', str(self.root/'missing')],
                                capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, b'')

    def test_helper_dependency_change_is_detected(self):
        source = self.stage / 'producer.C'
        header = self.stage / 'physics helper.h'
        source.write_text('#include "physics helper.h"\n')
        header.write_text('int constant=1;\n')
        depfile = self.stage / 'producer.d'
        escaped = str(header).replace(' ', '\\ ')
        depfile.write_text(f'producer: {source} ' + escaped + '\n')
        old = pub.dependency_identity([depfile])
        header.write_text('int constant=2;\n')
        self.assertNotEqual(pub.dependency_identity([depfile]), old)

    def test_cross_manifest_ownership_refuses_later_overwrite(self):
        pub.publish(self.stage, self.archive, self.manifest, self.pairs, 'run1')
        second_stage = self.root/'stage2'
        second_stage.mkdir()
        (second_stage/'0.root').write_bytes(b'run2')
        other_manifest = self.root/'other/smearing_latest.json'
        with self.assertRaisesRegex(RuntimeError, 'cross-manifest'):
            pub.publish(second_stage, self.root/'archives/run2', other_manifest,
                        [(second_stage/'0.root', self.pairs[0][1])], 'run2')
        self.assertEqual(self.pairs[0][1].read_bytes(), b'new-0')
        self.assertFalse(other_manifest.exists())
        self.assertTrue((self.root/'archives/run2/0.root').is_file())

    def test_protected_input_and_parent_aliases_are_rejected(self):
        source = self.root/'raw.root'
        source.write_bytes(b'raw input')
        alias = self.root/'alias'
        alias.symlink_to(self.output, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'protected'):
            pub.publish(self.stage, self.archive, self.manifest,
                        [(self.stage/'0.root', source)], 'run1', [source])
        self.assertEqual(source.read_bytes(), b'raw input')
        self.assertFalse(self.archive.exists())
        with self.assertRaisesRegex(ValueError, 'distinct'):
            pub.publish(self.stage, self.archive, self.manifest,
                        [(self.stage/'0.root', self.output/'0.root'),
                         (self.stage/'1.root', alias/'0.root')], 'run1')


if __name__ == '__main__':
    unittest.main()
