import os
import sys
import unittest
import shutil
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.io.io_to_formats import _write_tabular
from gwaslab.io.io_compress import PANDAS_GZIP_COMPRESSION
from gwaslab.info.g_Log import Log
from gwaslab.algorithm.core.conversions import mlog10p_to_p_format_str


def make_sumstats_indels():
    rows = []
    # SNP: NEA and EA length equal (1)
    rows.append({
        "CHR": 1,
        "POS": 100,
        "EA": "G",
        "NEA": "A",
        "SNPID": "1:100_A_G",
        "rsID": "rs1",
        "STATUS": "9960099",
        "EAF": 0.2,
        "BETA": 0.01,
        "SE": 0.02,
        "P": 1e-5,
        "N": 10000,
    })
    # Insertion: EA length > 1, NEA length == 1
    rows.append({
        "CHR": 1,
        "POS": 200,
        "EA": "ATC",
        "NEA": "A",
        "SNPID": "1:200_A_ATC",
        "rsID": "rs2",
        "STATUS": "9960099",
        "EAF": 0.3,
        "BETA": -0.02,
        "SE": 0.03,
        "P": 2e-6,
        "N": 8000,
    })
    # Deletion: EA length == 1, NEA length > 1
    rows.append({
        "CHR": 1,
        "POS": 300,
        "EA": "A",
        "NEA": "ATC",
        "SNPID": "1:300_ATC_A",
        "rsID": "rs3",
        "STATUS": "9960099",
        "EAF": 0.4,
        "BETA": 0.05,
        "SE": 0.04,
        "P": 3e-7,
        "N": 12000,
    })
    return pd.DataFrame(rows)


def make_p_export_sumstats(include_mlog10p=True):
    row_underflow = {
        "CHR": 1,
        "POS": 100,
        "EA": "A",
        "NEA": "G",
        "SNPID": "rs_underflow",
        "P": 0.0,
    }
    row_recoverable = {
        "CHR": 1,
        "POS": 200,
        "EA": "G",
        "NEA": "A",
        "SNPID": "rs_recoverable",
        "P": 0.0,
    }
    row_normal = {
        "CHR": 1,
        "POS": 300,
        "EA": "C",
        "NEA": "T",
        "SNPID": "rs_normal",
        "P": 1e-5,
    }
    if include_mlog10p:
        row_underflow["MLOG10P"] = 350.0
        row_recoverable["MLOG10P"] = 10.0
        row_normal["MLOG10P"] = 5.0
    return pd.DataFrame([row_underflow, row_recoverable, row_normal])


def read_p_values_from_tsv(path):
    df = pd.read_csv(path, sep="\t", dtype={"P": str})
    return dict(zip(df["SNPID"], df["P"]))


class TestIOToFormats(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="tmp_io_to_formats_")
        df = make_sumstats_indels()
        self.gl = Sumstats(sumstats=df, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", rsid="rsID", p="P", verbose=False)
        self.gl.set_build("19", verbose=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_to_format_bed_positions_and_fields(self):
        prefix = os.path.join(self.tmpdir, "out_bed")
        self.gl.to_format(prefix, fmt="bed", gzip=False, verbose=False)
        out_path = prefix + ".bed"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            lines = [line.strip().split("\t") for line in f.readlines()]
        # SNP row: start = POS-1, end = POS-1 + len(NEA)
        chrom, start, end, alleles, strand, snpid = lines[0]
        self.assertEqual(int(start), 99)
        self.assertEqual(int(end), 100)
        self.assertEqual(alleles, "A/G")
        self.assertEqual(strand, "+")
        self.assertEqual(snpid, "1:100_A_G")
        # Insertion row: start=end=POS, allele "-/<EA_without_first>"
        _, start_ins, end_ins, alleles_ins, _, _ = lines[1]
        self.assertEqual(int(start_ins), 200)
        self.assertEqual(int(end_ins), 200)
        self.assertEqual(alleles_ins, "-/TC")
        # Deletion row: end = POS + len(NEA) - 1, allele "<NEA_without_first>/-"
        _, start_del, end_del, alleles_del, _, _ = lines[2]
        self.assertEqual(int(start_del), 300)
        self.assertEqual(int(end_del), 302)
        self.assertEqual(alleles_del, "TC/-")

    def test_to_format_vep_positions_and_fields(self):
        prefix = os.path.join(self.tmpdir, "out_vep")
        self.gl.to_format(prefix, fmt="vep", gzip=False, verbose=False)
        out_path = prefix + ".vep"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            lines = [line.strip().split("\t") for line in f.readlines()]
        # SNP row: start=end=POS + len(NEA) - 1
        _, start, end, alleles, strand, snpid = lines[0]
        self.assertEqual(int(start), 100)
        self.assertEqual(int(end), 100)
        self.assertEqual(alleles, "A/G")
        self.assertEqual(strand, "+")
        self.assertEqual(snpid, "1:100_A_G")
        # Insertion row: start = POS+1, end = POS, allele "-/<EA_without_first>"
        _, start_ins, end_ins, alleles_ins, _, _ = lines[1]
        self.assertEqual(int(start_ins), 201)
        self.assertEqual(int(end_ins), 200)
        self.assertEqual(alleles_ins, "-/TC")
        # Deletion row: start = POS+1, end = POS + len(NEA) - 1, allele "<NEA_without_first>/-"
        _, start_del, end_del, alleles_del, _, _ = lines[2]
        self.assertEqual(int(start_del), 301)
        self.assertEqual(int(end_del), 302)
        self.assertEqual(alleles_del, "TC/-")

    def test_to_format_annovar_positions_and_fields(self):
        prefix = os.path.join(self.tmpdir, "out_annovar")
        self.gl.to_format(prefix, fmt="annovar", gzip=False, verbose=False)
        out_path = prefix + ".annovar"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            lines = [line.strip().split("\t") for line in f.readlines()]
        # SNP row: start = POS, end = POS-1 + len(NEA), columns NEA_out, EA_out
        _, start, end, nea_out, ea_out, snpid = lines[0]
        self.assertEqual(int(start), 100)
        self.assertEqual(int(end), 100)
        self.assertEqual(nea_out, "A")
        self.assertEqual(ea_out, "G")
        self.assertEqual(snpid, "1:100_A_G")
        # Insertion row: start=end=POS (ANNOVAR convention), NEA_out="-", EA_out without first base
        _, start_ins, end_ins, nea_out_ins, ea_out_ins, _ = lines[1]
        self.assertEqual(int(start_ins), 200)
        self.assertEqual(int(end_ins), 200)
        self.assertEqual(nea_out_ins, "-")
        self.assertEqual(ea_out_ins, "TC")
        # Deletion row: start=POS, end=POS-1 + len(NEA), NEA_out removed first base, EA_out="-"
        _, start_del, end_del, nea_out_del, ea_out_del, _ = lines[2]
        self.assertEqual(int(start_del), 300)
        self.assertEqual(int(end_del), 302)
        self.assertEqual(nea_out_del, "TC")
        self.assertEqual(ea_out_del, "-")

    def test_to_format_vcf_file_and_header(self):
        prefix = os.path.join(self.tmpdir, "out_vcf")
        self.gl.to_format(prefix, fmt="vcf", verbose=False)
        out_path = prefix + ".vcf"
        self.assertTrue(os.path.exists(out_path))
        with open(out_path, "r") as f:
            content = f.read()
        self.assertIn("##gwaslab_version", content)
        self.assertIn("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", content)

    def test_to_format_ssf_gzip_exists(self):
        prefix = os.path.join(self.tmpdir, "out_ssf")
        self.gl.to_format(prefix, fmt="ssf", tab_fmt="tsv", gzip=True, ssfmeta=True, verbose=False)
        out_path = prefix + ".ssf.tsv.gz"
        self.assertTrue(os.path.exists(out_path))
        self.assertGreater(os.path.getsize(out_path), 0)

    def test_write_tabular_gzip_uses_compresslevel_6(self):
        df = make_sumstats_indels()
        captured = {}

        def fake_to_csv(self, path, index=False, **kwargs):
            captured["path"] = path
            captured["kwargs"] = kwargs

        old_to_csv = pd.DataFrame.to_csv
        pd.DataFrame.to_csv = fake_to_csv
        try:
            to_csvargs = {}
            _write_tabular(
                df,
                {"CHR": "CHR"},
                os.path.join(self.tmpdir, "out.ssf.tsv.gz"),
                "tsv",
                to_csvargs,
                {},
                Log(),
                False,
                gzip=True,
            )
            self.assertEqual(
                captured["kwargs"]["compression"],
                PANDAS_GZIP_COMPRESSION,
            )
        finally:
            pd.DataFrame.to_csv = old_to_csv

    def test_write_tabular_respects_user_compression(self):
        df = make_sumstats_indels()
        captured = {}

        def fake_to_csv(self, path, index=False, **kwargs):
            captured["kwargs"] = kwargs

        old_to_csv = pd.DataFrame.to_csv
        pd.DataFrame.to_csv = fake_to_csv
        try:
            user_compression = {"method": "gzip", "compresslevel": 1}
            to_csvargs = {"compression": user_compression}
            _write_tabular(
                df,
                {"CHR": "CHR"},
                os.path.join(self.tmpdir, "out.ssf.tsv.gz"),
                "tsv",
                to_csvargs,
                {},
                Log(),
                False,
                gzip=True,
            )
            self.assertIs(captured["kwargs"]["compression"], user_compression)
        finally:
            pd.DataFrame.to_csv = old_to_csv


class TestLDAKFormat(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="tmp_io_ldak_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _export(self, name, statistics):
        data = {
            "SNPID": ["variant_1", "variant_2"],
            "CHR": [1, 2],
            "POS": [100, 200],
            "EA": ["G", "T"],
            "NEA": ["A", "C"],
            "N": [1000, 900],
        }
        data.update(statistics)
        column_args = {
            "snpid": "SNPID",
            "chrom": "CHR",
            "pos": "POS",
            "ea": "EA",
            "nea": "NEA",
            "n": "N",
        }
        for column, argument in {
            "Z": "z",
            "BETA": "beta",
            "SE": "se",
            "EAF": "eaf",
        }.items():
            if column in data:
                column_args[argument] = column
        sumstats = Sumstats(sumstats=pd.DataFrame(data), verbose=False, **column_args)
        prefix = os.path.join(self.tmpdir, name)
        sumstats.to_format(prefix, fmt="ldak", gzip=False, verbose=False)
        return prefix + ".ldak.tsv"

    def test_ldak_export_z_header_order_and_a1freq(self):
        out_path = self._export(
            "z",
            {
                "Z": [2.5, -1.5],
                "EAF": [0.25, 0.4],
            },
        )

        exported = pd.read_csv(out_path, sep="\t")
        self.assertListEqual(
            exported.columns.tolist(),
            ["Predictor", "A1", "A2", "Z", "n", "A1Freq"],
        )
        self.assertListEqual(exported["Predictor"].tolist(), ["variant_1", "variant_2"])

    def test_ldak_export_beta_se_header_order_without_a1freq(self):
        out_path = self._export(
            "beta_se",
            {
                "BETA": [0.1, -0.2],
                "SE": [0.01, 0.02],
            },
        )

        exported = pd.read_csv(out_path, sep="\t")
        self.assertListEqual(
            exported.columns.tolist(),
            ["Predictor", "A1", "A2", "BETA", "SE", "n"],
        )
        self.assertNotIn("A1Freq", exported.columns)

    def test_ldak_roundtrip_preserves_generic_predictor_as_snpid(self):
        out_path = self._export(
            "roundtrip",
            {
                "Z": [2.5, -1.5],
                "EAF": [0.25, 0.4],
            },
        )

        imported = Sumstats(sumstats=out_path, fmt="ldak", verbose=False)
        expected_columns = {"SNPID", "EA", "NEA", "Z", "N", "EAF"}
        self.assertTrue(expected_columns.issubset(imported.data.columns))
        self.assertNotIn("rsID", imported.data.columns)
        self.assertListEqual(
            imported.data["SNPID"].astype(str).tolist(),
            ["variant_1", "variant_2"],
        )
        self.assertListEqual(imported.data["N"].tolist(), [1000, 900])


class TestPExportUnderflow(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="tmp_io_p_underflow_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _export_p_values(self, df, **to_format_kwargs):
        gl = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            p="P",
            verbose=False,
        )
        prefix = os.path.join(self.tmpdir, "out_p")
        gl.to_format(prefix, fmt="gwaslab", gzip=False, verbose=False, **to_format_kwargs)
        return read_p_values_from_tsv(prefix + ".gwaslab.tsv")

    def test_p_underflow_export_from_mlog10p(self):
        p_values = self._export_p_values(make_p_export_sumstats())
        self.assertEqual(p_values["rs_underflow"], "1.0000e-350")
        self.assertNotEqual(p_values["rs_underflow"], "0.0000e+00")

    def test_p_recoverable_mlog10p_export(self):
        p_values = self._export_p_values(make_p_export_sumstats())
        self.assertEqual(p_values["rs_recoverable"], "1.0000e-10")

    def test_p_normal_export_unchanged(self):
        p_values = self._export_p_values(make_p_export_sumstats())
        self.assertEqual(p_values["rs_normal"], "1.0000e-05")

    def test_p_zero_without_mlog10p_stays_zero(self):
        p_values = self._export_p_values(make_p_export_sumstats(include_mlog10p=False))
        self.assertEqual(p_values["rs_underflow"], "0.0000e+00")
        self.assertEqual(p_values["rs_recoverable"], "0.0000e+00")

    def test_p_custom_float_format_from_mlog10p(self):
        df = pd.DataFrame({
            "CHR": [1],
            "POS": [100],
            "EA": ["A"],
            "NEA": ["G"],
            "SNPID": ["rs_custom"],
            "P": [0.0],
            "MLOG10P": [50.0],
        })
        p_values = self._export_p_values(
            df,
            float_formats={"P": "{:.2e}"},
        )
        self.assertEqual(p_values["rs_custom"], "1.00e-50")
        self.assertEqual(p_values["rs_custom"], mlog10p_to_p_format_str(50.0, "{:.2e}"))


if __name__ == "__main__":
    unittest.main()
