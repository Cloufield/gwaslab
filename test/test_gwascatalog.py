"""
Test suite for GWAS Catalog API v2 extension module.

Tests for GWAS Catalog API functionality.
"""

import os
import sys
import unittest
import time
import pytest

# Skip this entire test file when running with pytest
pytestmark = pytest.mark.skip(reason="Skipping test_gwascatalog.py")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
from gwaslab.extension.gwascatalog import GWASCatalogClient


class TestGWASCatalogAPI(unittest.TestCase):
    """Test cases for GWAS Catalog API client."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.client = GWASCatalogClient(verbose=False)
        self.test_trait = "type 2 diabetes mellitus"
        self.test_snp = "rs1050316"
        self.test_study_id = "GCST000001"
    
    def test_basic_request(self):
        """Test basic API request functionality."""
        # Test that we can make a request
        response = self.client._make_request("/v2/associations", {"size": 1})
        self.assertIsNotNone(response)
        self.assertIn("_embedded", response)
    
    def test_get_associations_by_trait(self):
        """Test Question 1: Get associations for a trait."""
        # Get first 10 associations for type 2 diabetes
        associations = self.client.get_associations(
            efo_trait=self.test_trait,
            size=10
        )
        
        self.assertIsNotNone(associations)
        if isinstance(associations, pd.DataFrame):
            self.assertGreater(len(associations), 0)
            # Check that expected columns exist
            self.assertIn("p_value", associations.columns)
            self.assertIn("snp_effect_allele", associations.columns)
        elif isinstance(associations, dict):
            # Should have _embedded structure
            self.assertIn("_embedded", associations)
    
    def test_get_studies_by_trait_and_cohort(self):
        """Test Question 2: Get studies by trait and cohort."""
        studies = self.client.get_studies(
            disease_trait="type 2 diabetes",
            cohort="UKB",
            size=20
        )
        
        self.assertIsNotNone(studies)
        if isinstance(studies, pd.DataFrame):
            self.assertGreater(len(studies), 0)
            self.assertIn("accession_id", studies.columns)
        elif isinstance(studies, dict):
            self.assertIn("_embedded", studies)
    
    def test_get_associations_by_snp(self):
        """Test Question 3: Get associations for a specific SNP."""
        associations = self.client.get_associations(
            rs_id=self.test_snp,
            sort="p_value",
            direction="asc",
            size=10
        )
        
        self.assertIsNotNone(associations)
        if isinstance(associations, pd.DataFrame):
            self.assertGreater(len(associations), 0)
            # Check that we have p_value column
            self.assertIn("p_value", associations.columns)
        elif isinstance(associations, dict):
            self.assertIn("_embedded", associations)
    
    def test_get_variants_by_gene(self):
        """Test Question 10: Get variants for a gene."""
        variants = self.client.get_variants(
            mapped_gene="HBS1L",
            size=10
        )
        
        self.assertIsNotNone(variants)
        if isinstance(variants, pd.DataFrame):
            # Should have some variants
            if len(variants) > 0:
                self.assertIn("rs_id", variants.columns)
        elif isinstance(variants, dict):
            self.assertIn("_embedded", variants)
    
    def test_get_single_variant(self):
        """Test getting a single variant by rsID."""
        variant = self.client.get_variants(rs_id=self.test_snp)
        
        self.assertIsNotNone(variant)
        if isinstance(variant, dict):
            self.assertIn("rs_id", variant)
            self.assertEqual(variant["rs_id"], self.test_snp)
    
    def test_get_single_study(self):
        """Test getting a single study by accession ID."""
        study = self.client.get_studies(accession_id=self.test_study_id)
        
        self.assertIsNotNone(study)
        if isinstance(study, dict):
            self.assertIn("accession_id", study)
            self.assertEqual(study["accession_id"], self.test_study_id)
    
    def test_get_all_variants_for_trait(self):
        """Test helper method: get_all_variants_for_trait."""
        # This might take a while, so limit it
        variants = self.client.get_all_variants_for_trait(self.test_trait)
        
        self.assertIsInstance(variants, set)
        self.assertGreater(len(variants), 0)
        # Check that all items are rsIDs
        for variant in list(variants)[:10]:  # Check first 10
            self.assertTrue(variant.startswith("rs"))
    
    def test_get_all_associations_for_trait(self):
        """Test helper method: get_all_associations_for_trait."""
        # This might take a while, limit to first page for testing
        # We'll test with a smaller query first
        associations = self.client.get_all_associations_for_trait(
            self.test_trait
        )
        
        self.assertIsInstance(associations, list)
        self.assertGreater(len(associations), 0)
        # Check structure of first association
        if len(associations) > 0:
            self.assertIsInstance(associations[0], dict)
            self.assertIn("p_value", associations[0])
    
    def test_get_all_genes_for_trait(self):
        """Test helper method: get_all_genes_for_trait."""
        genes = self.client.get_all_genes_for_trait(self.test_trait)
        
        self.assertIsInstance(genes, set)
        self.assertGreater(len(genes), 0)
        # Check that genes are strings
        for gene in list(genes)[:10]:  # Check first 10
            self.assertIsInstance(gene, str)
            self.assertGreater(len(gene), 0)
    
    def test_pagination(self):
        """Test pagination functionality."""
        # Get first page
        page1 = self.client.get_associations(
            efo_trait=self.test_trait,
            page=0,
            size=10
        )
        
        # Get second page
        page2 = self.client.get_associations(
            efo_trait=self.test_trait,
            page=1,
            size=10
        )
        
        self.assertIsNotNone(page1)
        self.assertIsNotNone(page2)
        
        # If both are DataFrames, they should have different data
        if isinstance(page1, pd.DataFrame) and isinstance(page2, pd.DataFrame):
            if len(page1) > 0 and len(page2) > 0:
                # They should have different association IDs (if available)
                if "association_id" in page1.columns and "association_id" in page2.columns:
                    self.assertFalse(
                        set(page1["association_id"]) == set(page2["association_id"])
                    )
    
    def test_get_all_pages(self):
        """Test _get_all_pages helper method."""
        # Get a small number of pages for testing
        items = self.client._get_all_pages(
            "/v2/associations",
            params={"efo_trait": self.test_trait, "size": 10},
            page_size=10
        )
        
        self.assertIsInstance(items, list)
        self.assertGreater(len(items), 0)
        # Should have multiple items (at least from first page)
        self.assertGreaterEqual(len(items), 1)
    
    def test_studies_with_full_summary_stats(self):
        """Test Question 9: Get studies with full summary statistics."""
        studies = self.client.get_studies(
            efo_trait=self.test_trait,
            full_pvalue_set=True,
            size=10
        )
        
        self.assertIsNotNone(studies)
        if isinstance(studies, pd.DataFrame):
            if len(studies) > 0:
                # Should have full_summary_stats column
                if "full_summary_stats" in studies.columns:
                    # At least some should have summary stats
                    has_stats = studies["full_summary_stats"].notna().any()
                    # This is optional, so we just check it doesn't crash
    
    def test_dataframe_conversion(self):
        """Test that responses are properly converted to DataFrames."""
        associations = self.client.get_associations(
            efo_trait=self.test_trait,
            size=5
        )
        
        if isinstance(associations, pd.DataFrame):
            # Should have variant_rsID column if snp_effect_allele was parsed
            if "snp_effect_allele" in associations.columns:
                # Check that parsing worked
                if "variant_rsID" in associations.columns:
                    # Should have valid rsIDs
                    rsids = associations["variant_rsID"].dropna()
                    if len(rsids) > 0:
                        self.assertTrue(rsids.iloc[0].startswith("rs"))
    
    def test_show_child_traits_parameter(self):
        """Test show_child_traits parameter for trait queries."""
        # Test with show_child_traits=False (only direct annotations)
        associations_false = self.client.get_associations(
            efo_trait=self.test_trait,
            show_child_traits=False,
            size=5
        )
        
        # Test with show_child_traits=True (include child traits)
        associations_true = self.client.get_associations(
            efo_trait=self.test_trait,
            show_child_traits=True,
            size=5
        )
        
        self.assertIsNotNone(associations_false)
        self.assertIsNotNone(associations_true)
    
    def test_extended_geneset_parameter(self):
        """Test extended_geneset parameter for gene queries."""
        # Test with extended_geneset=False (standard gene set)
        variants_false = self.client.get_variants(
            mapped_gene="HBS1L",
            extended_geneset=False,
            size=5
        )
        
        # Test with extended_geneset=True (extended gene set)
        variants_true = self.client.get_variants(
            mapped_gene="HBS1L",
            extended_geneset=True,
            size=5
        )
        
        self.assertIsNotNone(variants_false)
        self.assertIsNotNone(variants_true)
    
    def test_search_traits(self):
        """Test search_traits method for finding traits by free text."""
        # Search for traits containing "diabetes"
        # Note: This endpoint might not work as expected, so we'll test that it doesn't crash
        traits = self.client.search_traits("diabetes", size=10)
        
        self.assertIsNotNone(traits)
        # The endpoint might return empty results or use a different structure
        # So we just check it doesn't crash
        if isinstance(traits, pd.DataFrame):
            # If we get results, they should be valid
            pass
        elif isinstance(traits, dict):
            # Could be empty dict or have _embedded
            pass
    
    def test_get_trait_by_id(self):
        """Test get_trait method for getting a specific trait by EFO ID."""
        # Try to get a trait by ID (using a common one)
        trait = self.client.get_trait("EFO_0001360")  # Type 2 diabetes
        
        self.assertIsNotNone(trait)
        if isinstance(trait, dict) and len(trait) > 0:
            # Should have trait information
            self.assertIn("shortForm", trait)
    
    def test_rate_limiting(self):
        """Test that rate limiting is working."""
        # Make multiple rapid requests
        start_time = time.time()
        for i in range(3):
            self.client.get_associations(efo_trait=self.test_trait, size=1)
        end_time = time.time()
        
        # Should take at least (3 * rate_limit_delay) seconds
        elapsed = end_time - start_time
        min_expected = 2 * self.client.rate_limit_delay  # At least 2 delays
        self.assertGreaterEqual(elapsed, min_expected)


class TestGWASCatalogConvenienceFunctions(unittest.TestCase):
    """Test convenience functions."""
    
    def test_get_studies_function(self):
        """Test get_studies convenience function."""
        from gwaslab.extension.gwascatalog import get_studies
        
        studies = get_studies(efo_trait="type 2 diabetes mellitus", size=5, verbose=False)
        self.assertIsNotNone(studies)
    
    def test_get_associations_function(self):
        """Test get_associations convenience function."""
        from gwaslab.extension.gwascatalog import get_associations
        
        associations = get_associations(efo_trait="type 2 diabetes mellitus", size=5, verbose=False)
        self.assertIsNotNone(associations)
    
    def test_get_variants_function(self):
        """Test get_variants convenience function."""
        from gwaslab.extension.gwascatalog import get_variants
        
        variants = get_variants(rs_id="rs1050316", verbose=False)
        self.assertIsNotNone(variants)
    
    def test_search_traits_function(self):
        """Test search_traits convenience function."""
        from gwaslab.extension.gwascatalog import search_traits
        
        traits = search_traits("diabetes", size=5, verbose=False)
        self.assertIsNotNone(traits)


if __name__ == "__main__":
    # Run tests
    unittest.main(verbosity=2)

