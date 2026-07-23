"""Offline tests: get_novel forwards GWAS Catalog query params."""

import os
import sys
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.extension.gwascatalog.gwascatalog import (
    _merge_catalog_query_params,
    GWASCatalogClient,
)
from gwaslab.g_Sumstats import Sumstats


def _make_sumstats_one_lead():
    df = pd.DataFrame({
        "CHR": [1],
        "POS": [5000000],
        "P": [1e-10],
        "MLOG10P": [10.0],
        "SNPID": ["1:5000000_A_G"],
        "EA": ["A"],
        "NEA": ["G"],
    })
    ss = Sumstats(
        sumstats=df,
        chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
        snpid="SNPID", ea="EA", nea="NEA",
        verbose=False,
    )
    ss.build = "38"
    return ss


class TestMergeCatalogQueryParams(unittest.TestCase):
    def test_reserved_key_raises(self):
        with self.assertRaises(ValueError) as ctx:
            _merge_catalog_query_params(catalog_kwargs={"efo_id": "EFO_0006335"})
        self.assertIn("efo_id", str(ctx.exception))

    def test_explicit_sort_overlays_catalog_kwargs(self):
        merged = _merge_catalog_query_params(
            size=500,
            sort="p_value",
            direction="asc",
            catalog_kwargs={"extended_geneset": True, "sort": "or_value"},
        )
        self.assertEqual(merged["size"], 500)
        self.assertEqual(merged["sort"], "p_value")
        self.assertEqual(merged["direction"], "asc")
        self.assertTrue(merged["extended_geneset"])


class TestGetNovelGwascatalogPassthrough(unittest.TestCase):
    @patch("gwaslab.util.util_in_get_sig.GWASCatalogClient")
    def test_forwards_size_and_catalog_kwargs(self, mock_client_cls):
        mock_client = MagicMock()
        mock_client_cls.return_value = mock_client
        mock_client.get_known_variants_for_trait.return_value = pd.DataFrame({
            "SNPID": ["rs1"],
            "CHR": [1],
            "POS": [1000000],
        })

        ss = _make_sumstats_one_lead()
        ss.get_novel(
            efo="EFO_0006335",
            build="38",
            size=500,
            catalog_kwargs={"extended_geneset": True},
            verbose=False,
        )

        mock_client.get_known_variants_for_trait.assert_called_once()
        kwargs = mock_client.get_known_variants_for_trait.call_args.kwargs
        self.assertEqual(kwargs["size"], 500)
        self.assertIsNone(kwargs["sort"])
        self.assertEqual(kwargs["catalog_kwargs"], {"extended_geneset": True})

    @patch("gwaslab.util.util_in_get_sig.GWASCatalogClient")
    def test_forwards_sort_and_direction(self, mock_client_cls):
        mock_client = MagicMock()
        mock_client_cls.return_value = mock_client
        mock_client.get_known_variants_for_trait.return_value = pd.DataFrame({
            "SNPID": ["rs1"],
            "CHR": [1],
            "POS": [1000000],
        })

        ss = _make_sumstats_one_lead()
        ss.get_novel(
            efo="EFO_0006335",
            build="38",
            sort="p_value",
            direction="asc",
            verbose=False,
        )

        kwargs = mock_client.get_known_variants_for_trait.call_args.kwargs
        self.assertEqual(kwargs["sort"], "p_value")
        self.assertEqual(kwargs["direction"], "asc")

    def test_get_known_variants_for_trait_rejects_reserved_catalog_kwargs(self):
        client = GWASCatalogClient(verbose=False)
        with self.assertRaises(ValueError):
            client.get_known_variants_for_trait(
                efo="EFO_0006335",
                use_cache=False,
                catalog_kwargs={"page": 0},
                verbose=False,
            )


if __name__ == "__main__":
    unittest.main()
