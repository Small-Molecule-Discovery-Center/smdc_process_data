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
        report = Path(__file__).resolve().parent.parent / "data/masslynx_example_input.rpt"
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


if __name__ == "__main__":
    unittest.main()
