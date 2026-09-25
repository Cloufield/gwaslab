"""Issue #230: two VCF records at one indel site (not strand).

split_norm panels can list both REF=T ALT=TG AF=0 (insertion) and
REF=TG ALT=T AF=0.862 (deletion). Insertion is written first so current
check_unkonwn_indel first-match logic flips the reporter alleles.

Tests are written first (TDD). Production code should assign the sumstats
row to the VCF record with the smaller independent DAF.
"""
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
from pysam import VariantFile

from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm.hm_harmonize_sumstats import check_unkonwn_indel

CHROM = "1"
POS = 10
DELETION_AF = 0.862
INSERTION_AF = 0.0
KEEP_EAF = 0.845
KEEP_BETA = -0.035
SWAP_EAF = 0.155
SWAP_BETA = 0.035
DAF_TOLERANCE = 0.20
REF_MAF_THRESHOLD = 0.40


def _bgzip_tabix_vcf(vcf_path):
    bgzip_cmds = ["/home/yunye/tools/bin/bgzip", "bgzip"]
    tabix_cmds = ["/home/yunye/tools/bin/tabix", "tabix"]
    compressed = vcf_path + ".gz"
    bgzip_ok = False
    for cmd in bgzip_cmds:
        try:
            subprocess.run([cmd, "-f", vcf_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            bgzip_ok = True
            break
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    if not bgzip_ok:
        raise unittest.SkipTest("bgzip is required to index the two-record test VCF")
    tabix_ok = False
    for cmd in tabix_cmds:
        try:
            subprocess.run([cmd, "-p", "vcf", compressed], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            tabix_ok = True
            break
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    if not tabix_ok:
        raise unittest.SkipTest("tabix is required to index the two-record test VCF")
    return compressed


def _compress_fasta(fasta_path):
    bgzip_cmds = ["/home/yunye/tools/bin/bgzip", "bgzip"]
    for cmd in bgzip_cmds:
        try:
            subprocess.run([cmd, "-f", fasta_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return fasta_path + ".gz"
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    import gzip
    compressed = fasta_path + ".gz"
    with open(fasta_path, "rb") as src, gzip.open(compressed, "wb") as dst:
        dst.write(src.read())
    os.remove(fasta_path)
    return compressed


def write_two_record_vcf(directory, insertion_first=True):
    """Write REF=T/ALT=TG AF=0 and REF=TG/ALT=T AF=0.862. Insertion first by default."""
    vcf_path = os.path.join(directory, "two_records.vcf")
    insertion = (CHROM, POS, "ins", "T", "TG", INSERTION_AF)
    deletion = (CHROM, POS, "del", "TG", "T", DELETION_AF)
    rows = [insertion, deletion] if insertion_first else [deletion, insertion]
    with open(vcf_path, "w") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##contig=<ID=1,length=1000>\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, rsid, ref, alt, af in rows:
            handle.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\tAF={af}\n")
    return _bgzip_tabix_vcf(vcf_path)


def write_both_alleles_on_ref_fasta(directory):
    """FASTA so both T and TG match at POS (digit 6 = 6 after check_ref)."""
    fasta_path = os.path.join(directory, "ref.fa")
    seq = ["N"] * 100
    # 1-based POS=10 -> index 9: T then G so NEA=T and NEA=TG both match
    seq[POS - 1] = "T"
    seq[POS] = "G"
    for i, base in enumerate(seq):
        if base == "N":
            seq[i] = "A"
    with open(fasta_path, "w") as handle:
        handle.write(">1\n")
        handle.write("".join(seq) + "\n")
    return _compress_fasta(fasta_path)


def _digit7(status):
    return int(status) % 10


class TestIssue230TwoVcfRecords(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dir = self._tmpdir.name
        self.vcf_gz = write_two_record_vcf(self.dir, insertion_first=True)
        self.fasta_gz = write_both_alleles_on_ref_fasta(self.dir)

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_check_unkonwn_indel_keeps_direct_deletion(self):
        """Reporter alleles: NEA=TG EA=T EAF=0.845 must stay on the deletion record (digit 7=3)."""
        with VariantFile(self.vcf_gz) as vcf_reader:
            status = check_unkonwn_indel(
                CHROM,
                POS - 1,
                POS,
                "TG",
                "T",
                KEEP_EAF,
                vcf_reader,
                "AF",
                REF_MAF_THRESHOLD,
                1960369,
                None,
                DAF_TOLERANCE,
            )
        self.assertEqual(_digit7(status), 3, f"STATUS={status}: expected digit 7=3 (direct deletion), not first-match flip")

    def _harmonize(self, ea, nea, eaf, beta, sweep_mode):
        sumstats_path = os.path.join(self.dir, f"sumstats_{ea}_{nea}_{sweep_mode}.tsv")
        pd.DataFrame(
            {
                "SNPID": [f"{CHROM}:{POS}:{nea}:{ea}"],
                "CHR": [1],
                "POS": [POS],
                "EA": [ea],
                "NEA": [nea],
                "EAF": [eaf],
                "BETA": [beta],
                "SE": [0.0175],
                "P": [1e-5],
                "N": [275488],
            }
        ).to_csv(sumstats_path, sep="\t", index=False)
        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            beta="BETA",
            se="SE",
            n="N",
            build="19",
            verbose=False,
        )
        gl.basic_check(verbose=False)
        gl.harmonize(
            basic_check=False,
            ref_seq=self.fasta_gz,
            ref_infer=self.vcf_gz,
            ref_alt_freq="AF",
            threads=1,
            remove=False,
            verbose=False,
            sweep_mode=sweep_mode,
            infer_strand_kwargs={"daf_tolerance": DAF_TOLERANCE},
        )
        return gl.data.iloc[0]

    def test_harmonize_keep_deletion_nonsweep(self):
        row = self._harmonize("T", "TG", KEEP_EAF, KEEP_BETA, sweep_mode=False)
        self.assertEqual(row["EA"], "T")
        self.assertEqual(row["NEA"], "TG")
        self.assertAlmostEqual(float(row["EAF"]), KEEP_EAF, places=5)
        self.assertAlmostEqual(float(row["BETA"]), KEEP_BETA, places=5)
        self.assertEqual(_digit7(row["STATUS"]), 3)

    def test_harmonize_keep_deletion_sweep(self):
        row = self._harmonize("T", "TG", KEEP_EAF, KEEP_BETA, sweep_mode=True)
        self.assertEqual(row["EA"], "T")
        self.assertEqual(row["NEA"], "TG")
        self.assertAlmostEqual(float(row["EAF"]), KEEP_EAF, places=5)
        self.assertAlmostEqual(float(row["BETA"]), KEEP_BETA, places=5)
        self.assertEqual(_digit7(row["STATUS"]), 3)

    def test_harmonize_assign_to_deletion_from_insertion_coding_nonsweep(self):
        row = self._harmonize("TG", "T", SWAP_EAF, SWAP_BETA, sweep_mode=False)
        self.assertEqual(row["EA"], "T")
        self.assertEqual(row["NEA"], "TG")
        self.assertAlmostEqual(float(row["EAF"]), KEEP_EAF, places=5)
        self.assertAlmostEqual(float(row["BETA"]), KEEP_BETA, places=5)
        self.assertEqual(_digit7(row["STATUS"]), 4)

    def test_harmonize_assign_to_deletion_from_insertion_coding_sweep(self):
        row = self._harmonize("TG", "T", SWAP_EAF, SWAP_BETA, sweep_mode=True)
        self.assertEqual(row["EA"], "T")
        self.assertEqual(row["NEA"], "TG")
        self.assertAlmostEqual(float(row["EAF"]), KEEP_EAF, places=5)
        self.assertAlmostEqual(float(row["BETA"]), KEEP_BETA, places=5)
        self.assertEqual(_digit7(row["STATUS"]), 4)


if __name__ == "__main__":
    unittest.main()
