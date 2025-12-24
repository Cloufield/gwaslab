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


if __name__ == "__main__":
    unittest.main()
