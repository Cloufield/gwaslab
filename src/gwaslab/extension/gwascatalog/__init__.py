"""
GWAS Catalog API v2 Extension Module

This module provides an interface to interact with the GWAS Catalog REST API v2.
It allows querying studies, associations, variants, traits, and other GWAS Catalog data.

Example usage:
    from gwaslab.extension.gwascatalog import GWASCatalogClient
    
    client = GWASCatalogClient()
    studies = client.get_studies(efo_trait="EFO_0001360")
    associations = client.get_associations(study_id="GCST000001")
"""

from .gwascatalog import (
    GWASCatalogClient,
    get_studies,
    get_associations,
    get_variants,
    search_studies,
    search_associations,
    search_variants,
    search_traits,
)

__all__ = [
    "GWASCatalogClient",
    "get_studies",
    "get_associations",
    "get_variants",
    "search_studies",
    "search_associations",
    "search_variants",
    "search_traits",
]

