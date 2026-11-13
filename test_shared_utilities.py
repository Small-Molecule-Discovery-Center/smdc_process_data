import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import mysql_cnx
import raw_data_parsers


class MySQLConfigurationTests(unittest.TestCase):
    def test_local_and_remote_credentials_and_escaped_password(self):
        env = dict(mysql_user="test", mysql80_pwd="local@:/", mysql80_port="3307",
                   mysql_pwd="remote@:/", mysql_port="3306")
        with patch.dict(os.environ, env, clear=True), patch.object(mysql_cnx, "create_engine") as create:
            mysql_cnx.create_mysql_cnx()
            url = create.call_args.args[0]
            self.assertEqual((url.host, url.port, url.password, url.database),
                             ("localhost", 3307, "local@:/", "QB3_HiTS"))
            mysql_cnx.create_mysql_cnx(host="example.org", db="assay")
            url = create.call_args.args[0]
            self.assertEqual((url.host, url.port, url.password, url.database),
                             ("example.org", 3306, "remote@:/", "assay"))


class XevoTests(unittest.TestCase):
    def test_real_report_and_sample_archive(self):
        report = Path(__file__).parent / "data/masslynx_example_input.rpt"
        with tempfile.TemporaryDirectory() as d:
            normal = raw_data_parsers.extract_rpt_data(str(report), outdir=d)
            self.assertEqual(len(normal), 28)
            self.assertEqual(len(list((Path(d) / "samples").glob("*.txt"))), 28)
            self.assertEqual(normal.iloc[0]["Well"], "A,5")
            self.assertTrue((Path(d) / "samples/Sample_001_data.txt").read_bytes().startswith(b";Mass"))
        with tempfile.TemporaryDirectory() as d:
            samples = raw_data_parsers.extract_rpt_data(str(report), samples_only=True, outdir=d)
            self.assertEqual(list(samples.columns), ["Sample", "Plate", "Well", "Fname", "intensitymax", "massatmax"])
            self.assertEqual(len(samples), 28)
            self.assertTrue((Path(d) / "samples.zip").is_file())
            self.assertTrue((Path(d) / "samples/extracted_samples.csv").is_file())

    def _two_spectrum_report(self, d, first_procdesc):
        """Copy the example report, giving sample 1 an extra spectrum before the real one,
        the way newer Xevo reports do."""
        lines = (Path(__file__).parent / "data/masslynx_example_input.rpt").read_bytes().splitlines(keepends=True)
        start = next(i for i, l in enumerate(lines) if l.startswith(b"[SPECTRUM]"))
        mass = next(i for i in range(start, len(lines)) if lines[i].startswith(b";Mass"))
        end = next(i for i in range(mass, len(lines)) if lines[i].lstrip().startswith(b"}"))
        extra = []
        for l in lines[start:end + 2]:  # through the closing braces of [MS] and [SPECTRUM]
            if l.startswith(b"ProcDesc\t"):
                l = b"ProcDesc\t" + first_procdesc + b"\r\n"
            elif l.startswith(b"BPI\t"):
                l = b"BPI\t1\r\n"
            elif l.startswith(b"BPM\t"):
                l = b"BPM\t1.0000\r\n"
            extra.append(l)
        report = Path(d) / "two_spectra.rpt"
        report.write_bytes(b"".join(lines[:start] + extra + lines[start:]))
        return report, len(extra)

    def test_two_spectra_per_sample_uses_maxent(self):
        with tempfile.TemporaryDirectory() as d:
            report, _ = self._two_spectrum_report(d, b"Combine (81:93)")
            df = raw_data_parsers.extract_rpt_data(str(report), outdir=d)
            self.assertEqual(len(df), 28)
            self.assertEqual(df.iloc[0]["intensitymax"], 11036179)
            self.assertEqual(df.iloc[0]["massatmax"], 25247.0)

    def test_two_maxent_spectra_uses_the_last(self):
        with tempfile.TemporaryDirectory() as d:
            report, n_extra = self._two_spectrum_report(d, b"MaxEnt 1 copy")
            df = raw_data_parsers.extract_rpt_data(str(report), outdir=d)
            self.assertEqual(df.iloc[0]["intensitymax"], 11036179)
            with tempfile.TemporaryDirectory() as d2:
                plain = raw_data_parsers.extract_rpt_data(
                    str(Path(__file__).parent / "data/masslynx_example_input.rpt"), outdir=d2)
                self.assertEqual(df.iloc[0]["linenum_start"], plain.iloc[0]["linenum_start"] + n_extra)
                self.assertEqual((Path(d) / "samples/Sample_001_data.txt").read_bytes(),
                                 (Path(d2) / "samples/Sample_001_data.txt").read_bytes())


if __name__ == "__main__":
    unittest.main()
