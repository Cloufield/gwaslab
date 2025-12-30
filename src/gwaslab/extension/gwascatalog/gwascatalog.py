"""
GWAS Catalog REST API v2 Client

This module provides a comprehensive interface to interact with the GWAS Catalog REST API v2.
It supports querying studies, associations, variants, traits, genes, and other entities.

API Documentation: https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference
Based on sample notebook examples from the GWAS Catalog API documentation.
"""

import requests
import pandas as pd
import json
import time
from typing import Optional, Dict, List, Any, Union
from gwaslab.info.g_Log import Log


class GWASCatalogClient:
    """
    Client for interacting with the GWAS Catalog REST API v2.
    
    This class provides methods to query various GWAS Catalog endpoints including
    studies, associations, variants, traits, genes, and more.
    
    Attributes
    ----------
    base_url : str
        Base URL for the GWAS Catalog API
    session : requests.Session
        HTTP session for making requests
    timeout : int
        Request timeout in seconds
    log : Log
        Logging object for tracking operations
        
    Examples
    --------
    >>> client = GWASCatalogClient()
    >>> studies = client.get_studies(efo_trait="type 2 diabetes mellitus")
    >>> associations = client.get_associations(efo_trait="type 2 diabetes mellitus")
    """
    
    BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
    MAX_REQUESTS_PER_SECOND = 15  # API rate limit
    
    def __init__(self, timeout: int = 30, log: Optional[Log] = None, verbose: bool = True, 
                 rate_limit_delay: float = 0.067):
        """
        Initialize the GWAS Catalog API client.
        
        Parameters
        ----------
        timeout : int, optional
            Request timeout in seconds (default: 30)
        log : Log, optional
            Logging object for tracking operations (default: creates new Log)
        verbose : bool, optional
            Whether to print log messages (default: True)
        rate_limit_delay : float, optional
            Minimum delay between requests in seconds to respect rate limit (default: 0.067, ~15 req/s)
        """
        self.base_url = self.BASE_URL
        self.session = requests.Session()
        self.timeout = timeout
        self.verbose = verbose
        self.log = log if log is not None else Log()
        self.rate_limit_delay = rate_limit_delay
        self._last_request_time = 0.0
        
    def _make_request(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        suppress_request_logging: bool = False
    ) -> Optional[Dict[str, Any]]:
        """
        Make a GET request to the API endpoint with rate limiting.
        
        Parameters
        ----------
        endpoint : str
            API endpoint path (e.g., "/v2/studies", "/v2/associations")
        params : dict, optional
            Query parameters
        suppress_request_logging : bool, optional
            If True, suppress detailed request/response logging (default: False)
            
        Returns
        -------
        dict or None
            JSON response as a dictionary, or None if the request fails
        """
        # Rate limiting: ensure minimum delay between requests
        current_time = time.time()
        time_since_last = current_time - self._last_request_time
        if time_since_last < self.rate_limit_delay:
            sleep_time = self.rate_limit_delay - time_since_last
            time.sleep(sleep_time)
        self._last_request_time = time.time()
        
        full_url = f"{self.base_url}{endpoint}"
        
        try:
            if not suppress_request_logging:
                self.log.write(f"Requesting: {full_url}", verbose=self.verbose)
                if params:
                    self.log.write(f"Parameters: {params}", verbose=self.verbose)
                
            response = self.session.get(full_url, params=params, timeout=self.timeout)
            
            if response.status_code == 200:
                if not suppress_request_logging:
                    self.log.write(f"Status code: {response.status_code}", verbose=self.verbose)
                return response.json()
            else:
                # Always log errors
                self.log.warning(f"Error: Received status code {response.status_code}", verbose=self.verbose)
                self.log.warning(f"URL: {response.url}", verbose=self.verbose)
                self.log.warning(f"Response: {response.text}", verbose=self.verbose)
                return None
                
        except requests.exceptions.RequestException as e:
            # Always log errors
            self.log.warning(f"An error occurred: {e}", verbose=self.verbose)
            return None
    
    def _get_all_pages(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        page_size: int = 200,
        delay: float = 0.1
    ) -> List[Dict[str, Any]]:
        """
        Retrieve all pages of results from a paginated endpoint.
        
        Parameters
        ----------
        endpoint : str
            API endpoint path
        params : dict, optional
            Query parameters
        page_size : int, optional
            Number of items per page (default: 200)
        delay : float, optional
            Delay between requests in seconds (default: 0.1)
            
        Returns
        -------
        list
            List of all items from all pages
        """
        all_items = []
        current_page = 0
        total_pages = 1  # Initialize to 1 to start the loop
        
        # Set page size in params
        if params is None:
            params = {}
        params["size"] = page_size
        
        while current_page < total_pages:
            params["page"] = current_page
            
            # Suppress detailed request logging for pages 2 onwards
            suppress_logging = (current_page > 0)
            data = self._make_request(endpoint, params, suppress_request_logging=suppress_logging)
            
            if data and '_embedded' in data:
                # Extract items from embedded resources
                embedded_key = list(data['_embedded'].keys())[0] if data['_embedded'] else None
                if embedded_key:
                    items = data['_embedded'][embedded_key]
                    all_items.extend(items)
                
                # Update total pages from first response
                if current_page == 0:
                    total_pages = data.get('page', {}).get('totalPages', 1)
                    total_elements = data.get('page', {}).get('totalElements', 0)
                    self.log.write(f"Found {total_elements} total items across {total_pages} pages", verbose=self.verbose)
                
                # For page 1, show full details. For subsequent pages, show simplified progress
                if current_page == 0:
                    self.log.write(f"Page {current_page + 1}/{total_pages} processed. Collected {len(all_items)} items so far.", verbose=self.verbose)
                else:
                    # Show progress every 20 pages or on the last page, with simplified message
                    if (current_page + 1) % 20 == 0 or (current_page + 1) == total_pages:
                        self.log.write(f"Page {current_page + 1}/{total_pages}... ({len(all_items)} items)", verbose=self.verbose)
                
                current_page += 1
                
                if delay > 0:
                    time.sleep(delay)  # Be polite to the API
            else:
                self.log.warning("No more data or an error occurred. Stopping.", verbose=self.verbose)
                break
        
        self.log.write(f"Retrieved {len(all_items)} total items", verbose=self.verbose)
        return all_items
    
    # ============================================================================
    # Studies Endpoints
    # ============================================================================
    
    def get_studies(
        self,
        accession_id: Optional[str] = None,
        efo_trait: Optional[str] = None,
        disease_trait: Optional[str] = None,
        cohort: Optional[str] = None,
        full_pvalue_set: Optional[bool] = None,
        gxe: Optional[bool] = None,
        show_child_traits: Optional[bool] = None,
        page: Optional[int] = None,
        size: Optional[int] = None,
        get_all: bool = False
    ) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
        """
        Retrieve GWAS studies from the catalog.
        
        Parameters
        ----------
        accession_id : str, optional
            Specific study accession ID (e.g., "GCST000001")
        efo_trait : str, optional
            EFO trait name to filter studies
        disease_trait : str, optional
            Disease trait name to filter studies
        cohort : str, optional
            Cohort name to filter studies (e.g., "UKB")
        full_pvalue_set : bool, optional
            Filter for studies with full summary statistics
        gxe : bool, optional
            Filter for Gene-by-Environment (GxE) studies
        show_child_traits : bool, optional
            If True, include child traits in results (default: True per API).
            If False, only return data annotated directly with the query term.
        page : int, optional
            Page number (for pagination)
        size : int, optional
            Page size (for pagination)
        get_all : bool, optional
            If True, retrieve all pages (default: False)
            
        Returns
        -------
        pandas.DataFrame, dict, or list
            Study data as DataFrame (if multiple) or dict/list (if single)
        """
        if accession_id:
            endpoint = f"/v2/studies/{accession_id}"
            response = self._make_request(endpoint)
            return response if response else {}
        else:
            endpoint = "/v2/studies"
            params = {}
            if efo_trait:
                params["efo_trait"] = efo_trait
            if disease_trait:
                params["disease_trait"] = disease_trait
            if cohort:
                params["cohort"] = cohort
            if full_pvalue_set is not None:
                params["full_pvalue_set"] = full_pvalue_set
            if gxe is not None:
                params["gxe"] = gxe
            if show_child_traits is not None:
                params["show_child_traits"] = show_child_traits
                
            if get_all:
                items = self._get_all_pages(endpoint, params=params)
                return self._studies_to_dataframe(items) if items else pd.DataFrame()
            else:
                if page is not None:
                    params["page"] = page
                if size is not None:
                    params["size"] = size
                    
                response = self._make_request(endpoint, params=params)
                if response and '_embedded' in response:
                    items = response['_embedded']['studies']
                    return self._studies_to_dataframe(items) if items else pd.DataFrame()
                return response if response else {}
    
    # ============================================================================
    # Associations Endpoints
    # ============================================================================
    
    def get_associations(
        self,
        efo_trait: Optional[str] = None,
        rs_id: Optional[str] = None,
        accession_id: Optional[str] = None,
        sort: Optional[str] = None,
        direction: Optional[str] = None,
        show_child_traits: Optional[bool] = None,
        extended_geneset: Optional[bool] = None,
        page: Optional[int] = None,
        size: Optional[int] = None,
        get_all: bool = False
    ) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
        """
        Retrieve GWAS associations from the catalog.
        
        Parameters
        ----------
        efo_trait : str, optional
            EFO trait name to filter associations
        rs_id : str, optional
            Variant rsID to filter associations (e.g., "rs1050316")
        accession_id : str, optional
            Study accession ID to filter associations
        sort : str, optional
            Field to sort by (e.g., "p_value", "or_value")
        direction : str, optional
            Sort direction: "asc" or "desc"
        show_child_traits : bool, optional
            If True, include child traits in results (default: True per API).
            If False, only return data annotated directly with the query term.
        extended_geneset : bool, optional
            If True, use extended gene set (all Ensembl/RefSeq genes within 50kb).
            If False, use standard gene set (default: False).
        page : int, optional
            Page number (for pagination)
        size : int, optional
            Page size (for pagination)
        get_all : bool, optional
            If True, retrieve all pages (default: False)
            
        Returns
        -------
        pandas.DataFrame, dict, or list
            Association data as DataFrame (if multiple) or dict/list (if single)
        """
        endpoint = "/v2/associations"
        params = {}
        
        if efo_trait:
            params["efo_trait"] = efo_trait
        if rs_id:
            params["rs_id"] = rs_id
        if accession_id:
            params["accession_id"] = accession_id
        if sort:
            params["sort"] = sort
        if direction:
            params["direction"] = direction
        if show_child_traits is not None:
            params["show_child_traits"] = show_child_traits
        if extended_geneset is not None:
            params["extended_geneset"] = extended_geneset
            
        if get_all:
            items = self._get_all_pages(endpoint, params=params)
            return self._associations_to_dataframe(items) if items else pd.DataFrame()
        else:
            if page is not None:
                params["page"] = page
            if size is not None:
                params["size"] = size
                
            response = self._make_request(endpoint, params=params)
            if response and '_embedded' in response:
                items = response['_embedded']['associations']
                return self._associations_to_dataframe(items) if items else pd.DataFrame()
            return response if response else {}
    
    # ============================================================================
    # Variants (SNPs) Endpoints
    # ============================================================================
    
    def get_variants(
        self,
        rs_id: Optional[str] = None,
        mapped_gene: Optional[str] = None,
        extended_geneset: Optional[bool] = None,
        page: Optional[int] = None,
        size: Optional[int] = None,
        get_all: bool = False
    ) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
        """
        Retrieve genetic variants (SNPs) from the catalog.
        
        Parameters
        ----------
        rs_id : str, optional
            Specific variant rsID (e.g., "rs1050316")
        mapped_gene : str, optional
            Gene name to filter variants
        extended_geneset : bool, optional
            If True, use extended gene set (all Ensembl/RefSeq genes within 50kb).
            If False, use standard gene set (default: False).
        page : int, optional
            Page number (for pagination)
        size : int, optional
            Page size (for pagination)
        get_all : bool, optional
            If True, retrieve all pages (default: False)
            
        Returns
        -------
        pandas.DataFrame, dict, or list
            Variant data as DataFrame (if multiple) or dict/list (if single)
        """
        if rs_id:
            endpoint = f"/v2/single-nucleotide-polymorphisms/{rs_id}"
            response = self._make_request(endpoint)
            return response if response else {}
        else:
            endpoint = "/v2/single-nucleotide-polymorphisms"
            params = {}
            if mapped_gene:
                params["mapped_gene"] = mapped_gene
            if extended_geneset is not None:
                params["extended_geneset"] = extended_geneset
                
            if get_all:
                items = self._get_all_pages(endpoint, params=params)
                return self._variants_to_dataframe(items) if items else pd.DataFrame()
            else:
                if page is not None:
                    params["page"] = page
                if size is not None:
                    params["size"] = size
                    
                response = self._make_request(endpoint, params=params)
                if response and '_embedded' in response:
                    items = response['_embedded']['snps']
                    return self._variants_to_dataframe(items) if items else pd.DataFrame()
                return response if response else {}
    
    # ============================================================================
    # EFO Traits Endpoint
    # ============================================================================
    
    def search_traits(
        self,
        query: str,
        page: Optional[int] = None,
        size: Optional[int] = None,
        get_all: bool = False
    ) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
        """
        Search for EFO traits by free text query.
        
        This endpoint allows you to search for traits by name. For example, searching
        for "COVID-19" could return "COVID-19", "COVID-19 symptoms measurement",
        "long COVID-19", etc.
        
        Parameters
        ----------
        query : str
            Free text search term for trait name
        page : int, optional
            Page number (for pagination)
        size : int, optional
            Page size (for pagination)
        get_all : bool, optional
            If True, retrieve all pages (default: False)
            
        Returns
        -------
        pandas.DataFrame, dict, or list
            Trait data as DataFrame (if multiple) or dict/list (if single)
            
        Examples
        --------
        >>> client = GWASCatalogClient()
        >>> traits = client.search_traits("COVID-19")
        >>> traits = client.search_traits("diabetes", get_all=True)
        """
        endpoint = "/v2/efoTraits"
        params = {"query": query}
        
        if get_all:
            items = self._get_all_pages(endpoint, params=params)
            return self._traits_to_dataframe(items) if items else pd.DataFrame()
        else:
            if page is not None:
                params["page"] = page
            if size is not None:
                params["size"] = size
                
            response = self._make_request(endpoint, params=params)
            if response and '_embedded' in response:
                items = response['_embedded']['efoTraits']
                return self._traits_to_dataframe(items) if items else pd.DataFrame()
            return response if response else {}
    
    def get_trait(self, efo_id: str) -> Dict[str, Any]:
        """
        Get a specific EFO trait by its ID.
        
        Parameters
        ----------
        efo_id : str
            EFO trait ID (e.g., "MONDO_0004979", "EFO_0001360")
            
        Returns
        -------
        dict
            Trait information dictionary
        """
        endpoint = f"/v2/efoTraits/{efo_id}"
        response = self._make_request(endpoint)
        return response if response else {}
    
    # ============================================================================
    # Helper Methods for Common Queries
    # ============================================================================
    
    def get_all_variants_for_trait(self, trait_name: str) -> set:
        """
        Fetches all unique variant rsIDs for a given trait by handling API pagination.
        
        Parameters
        ----------
        trait_name : str
            The name of the trait to query
            
        Returns
        -------
        set
            A set of unique rsID strings associated with the trait
        """
        variants = set()
        self.log.write(f"--- Starting search for '{trait_name}' ---", verbose=self.verbose)
        
        associations = self.get_associations(efo_trait=trait_name, get_all=True)
        
        if isinstance(associations, pd.DataFrame):
            if 'snp_effect_allele' in associations.columns:
                for allele_list in associations['snp_effect_allele']:
                    if isinstance(allele_list, list) and len(allele_list) > 0:
                        risk_allele_str = allele_list[0]
                        rs_id = risk_allele_str.split('-')[0]
                        if rs_id.startswith('rs'):
                            variants.add(rs_id)
        elif isinstance(associations, list):
            for association in associations:
                if 'snp_effect_allele' in association and association['snp_effect_allele']:
                    risk_allele_str = association['snp_effect_allele'][0]
                    rs_id = risk_allele_str.split('-')[0]
                    if rs_id.startswith('rs'):
                        variants.add(rs_id)
        
        self.log.write(f"--- Finished fetching for '{trait_name}'. Found {len(variants)} total unique variants. ---", verbose=self.verbose)
        return variants
    
    def get_all_associations_for_trait(self, trait_name: str, sort_by_pvalue: bool = True) -> List[Dict[str, Any]]:
        """
        Fetches all association objects for a given trait by handling API pagination.
        
        Parameters
        ----------
        trait_name : str
            The name of the trait to query
        sort_by_pvalue : bool, optional
            If True, sort by p-value (default: True)
            
        Returns
        -------
        list
            A list of all association dictionaries for the trait
        """
        self.log.write(f"--- Starting search for all associations related to '{trait_name}' ---", verbose=self.verbose)
        
        params = {}
        if sort_by_pvalue:
            params["sort"] = "p_value"
            params["direction"] = "asc"
        
        items = self._get_all_pages("/v2/associations", params={**params, "efo_trait": trait_name}, page_size=40)
        
        self.log.write(f"--- Finished fetching. Found {len(items)} total associations. ---", verbose=self.verbose)
        return items
    
    def get_all_genes_for_trait(self, trait_name: str) -> set:
        """
        Fetches all unique gene names for a given trait.
        
        Parameters
        ----------
        trait_name : str
            The name of the trait to query
            
        Returns
        -------
        set
            A set of unique gene names associated with the trait
        """
        all_genes = set()
        associations = self.get_associations(efo_trait=trait_name, get_all=True)
        
        if isinstance(associations, pd.DataFrame):
            if 'mapped_genes' in associations.columns:
                for gene_list in associations['mapped_genes']:
                    if isinstance(gene_list, list):
                        for gene_string in gene_list:
                            genes = [gene.strip() for gene in gene_string.split(',') if gene.strip()]
                            all_genes.update(genes)
        elif isinstance(associations, list):
            for association in associations:
                if 'mapped_genes' in association and association['mapped_genes']:
                    for gene_string in association['mapped_genes']:
                        genes = [gene.strip() for gene in gene_string.split(',') if gene.strip()]
                        all_genes.update(genes)
        
        return all_genes
    
    # ============================================================================
    # Data Conversion Helpers
    # ============================================================================
    
    def _studies_to_dataframe(self, studies: List[Dict[str, Any]]) -> pd.DataFrame:
        """Convert studies list to pandas DataFrame."""
        if not studies:
            return pd.DataFrame()
        
        # Convert directly to DataFrame - the API returns flat structures
        df = pd.DataFrame(studies)
        return df
    
    def _associations_to_dataframe(self, associations: List[Dict[str, Any]]) -> pd.DataFrame:
        """Convert associations list to pandas DataFrame."""
        if not associations:
            return pd.DataFrame()
        
        # Convert directly to DataFrame - the API returns flat structures
        df = pd.DataFrame(associations)
        
        # Parse snp_effect_allele if present
        if 'snp_effect_allele' in df.columns:
            risk_allele_str = df['snp_effect_allele'].str[0]
            split_allele = risk_allele_str.str.split('-', n=1, expand=True)
            if len(split_allele.columns) >= 2:
                df['variant_rsID'] = split_allele[0]
                df['risk_allele_base'] = split_allele[1]
        
        return df
    
    def _variants_to_dataframe(self, variants: List[Dict[str, Any]]) -> pd.DataFrame:
        """Convert variants list to pandas DataFrame."""
        if not variants:
            return pd.DataFrame()
        
        # Convert directly to DataFrame
        df = pd.DataFrame(variants)
        return df
    
    def _traits_to_dataframe(self, traits: List[Dict[str, Any]]) -> pd.DataFrame:
        """Convert traits list to pandas DataFrame."""
        if not traits:
            return pd.DataFrame()
        
        # Convert directly to DataFrame
        df = pd.DataFrame(traits)
        return df
    
    def _extract_chr_pos_from_locations(self, locations: Any) -> tuple:
        """
        Extract chromosome and position from locations field.
        
        Parameters
        ----------
        locations : Any
            Locations value (usually a list like ['12:111803962'] or a string like "12:111803962")
            
        Returns
        -------
        tuple
            (chromosome, position) or (None, None) if extraction fails
        """
        if locations is None:
            return None, None
        
        # Convert to string if needed
        location_str = None
        
        # Handle if it's already a string (e.g., "12:111803962")
        if isinstance(locations, str):
            location_str = locations.strip()
        # Handle if it's a list
        elif isinstance(locations, list):
            if len(locations) == 0:
                return None, None
            # Get first element
            first_elem = locations[0]
            if isinstance(first_elem, str):
                location_str = first_elem.strip()
            elif isinstance(first_elem, (int, float)):
                # If it's just a number, we can't extract chr:pos
                return None, None
            else:
                # Try to convert to string
                location_str = str(first_elem).strip()
        # Handle other types - try to convert to string
        else:
            location_str = str(locations).strip()
        
        if not location_str or location_str == 'nan' or location_str.lower() == 'none':
            return None, None
        
        # Parse the location string
        # Format: "12:111803962" or "chr12:111803962" or "12:111803962-111803962" (range)
        # Remove 'chr' prefix if present
        location_str = location_str.replace('chr', '').replace('CHR', '')
        
        # Split by colon
        parts = location_str.split(':')
        if len(parts) >= 2:
            try:
                chr_val = parts[0].strip()
                # Handle position (might be a range like "111803962-111803963")
                pos_str = parts[1].strip()
                if '-' in pos_str:
                    # Take the first position if it's a range
                    pos_str = pos_str.split('-')[0].strip()
                pos_val = int(pos_str)
                
                # Validate chromosome (should be 1-22, X, Y, MT, or similar)
                if chr_val and pos_val > 0:
                    return chr_val, pos_val
            except (ValueError, IndexError):
                pass
        
        return None, None
    
    def get_known_variants_for_trait(
        self,
        efo: str,
        sig_level: float = 5e-8,
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Retrieve known variants from GWAS Catalog API v2 for a given EFO trait.
        
        This method handles MONDO to EFO conversion, retrieves associations,
        and returns them in the format expected by _get_novel.
        
        Parameters
        ----------
        efo : str
            EFO trait ID (e.g., "EFO_0001360"), MONDO ID (e.g., "MONDO_0005148"),
            or trait name (e.g., "duodenal ulcer")
        sig_level : float
            P-value threshold for filtering associations (default: 5e-8)
        verbose : bool
            Whether to print log messages (default: True)
            
        Returns
        -------
        pd.DataFrame
            DataFrame with columns: SNPID, CHR, POS, P, BETA, SE, OR, TRAIT,
            REPORT_GENENAME, PUBMEDID, AUTHOR, STUDY
        """
        self.log.write(" -Querying GWAS Catalog API v2 for trait: {}...".format(efo), verbose=verbose)
        
        try:
            # For MONDO IDs, first get the trait info to find EFO ID
            efo_to_use = efo
            trait_name_to_use = None
            
            if efo.startswith("MONDO_"):
                self.log.write(" -MONDO ID detected, looking up EFO equivalent...", verbose=verbose)
                try:
                    trait_info = self.get_trait(efo)
                    if trait_info:
                        if verbose:
                            self.log.write(" -Trait info keys: {}".format(list(trait_info.keys())[:10]), verbose=verbose)
                        
                        # Try to get EFO ID from various possible fields
                        efo_id = None
                        if 'shortForm' in trait_info and trait_info['shortForm']:
                            efo_id = trait_info['shortForm']
                        elif 'id' in trait_info and trait_info['id']:
                            if trait_info['id'].startswith('EFO_'):
                                efo_id = trait_info['id']
                            elif 'EFO_' in str(trait_info['id']):
                                # Extract EFO from ID string
                                parts = str(trait_info['id']).split('EFO_')
                                if len(parts) > 1:
                                    efo_id = 'EFO_' + parts[1].split('/')[0].split('_')[0]
                        
                        # Also try URI
                        if not efo_id and 'uri' in trait_info and trait_info['uri']:
                            uri = str(trait_info['uri'])
                            if 'EFO_' in uri:
                                parts = uri.split('EFO_')
                                if len(parts) > 1:
                                    efo_id = 'EFO_' + parts[1].split('/')[0].split('_')[0]
                        
                        # Get trait name as fallback
                        if 'trait' in trait_info:
                            trait_name_to_use = trait_info['trait']
                        
                        if efo_id and efo_id.startswith('EFO_'):
                            efo_to_use = efo_id
                            self.log.write(" -Found EFO ID: {} for MONDO ID: {}".format(efo_id, efo), verbose=verbose)
                        elif trait_name_to_use:
                            self.log.write(" -Using trait name '{}' instead of ID for MONDO: {}".format(trait_name_to_use, efo), verbose=verbose)
                            efo_to_use = trait_name_to_use
                        else:
                            self.log.write(" -Warning: Could not extract EFO ID or trait name. Trying MONDO ID directly...", verbose=verbose)
                    else:
                        self.log.write(" -Warning: Could not retrieve trait info for MONDO ID. Trying directly...", verbose=verbose)
                except Exception as e:
                    self.log.write(" -Warning: Error looking up trait info: {}. Trying MONDO ID directly...".format(str(e)), verbose=verbose)
            
            # Get associations from API v2 using the EFO ID or trait name
            associations_df = self.get_associations(efo_trait=efo_to_use, get_all=True, sort="p_value", direction="asc")
            
            # If no results and we tried a converted EFO ID, try the original MONDO ID with old API
            if (not isinstance(associations_df, pd.DataFrame) or len(associations_df) == 0) and efo.startswith("MONDO_"):
                self.log.write(" -No results from API v2, falling back to old API v1 for MONDO ID...", verbose=verbose)
                # Fallback to old API which supports MONDO IDs directly
                try:
                    from gwaslab.util.util_ex_gwascatalog import gwascatalog_trait
                    known_Sumstats = gwascatalog_trait(efo, source="NCBI", sig_level=sig_level, 
                                                       use_cache=True, cache_dir="./", 
                                                       verbose=verbose, log=self.log)
                    if known_Sumstats and hasattr(known_Sumstats, 'data') and len(known_Sumstats.data) > 0:
                        self.log.write(" -Retrieved {} variants using old API v1".format(len(known_Sumstats.data)), verbose=verbose)
                        return known_Sumstats.data.copy()
                except Exception as e2:
                    self.log.write(" -Fallback to old API also failed: {}".format(str(e2)), verbose=verbose)
            
            if not isinstance(associations_df, pd.DataFrame) or len(associations_df) == 0:
                self.log.write(" -No associations found for trait: {}".format(efo), verbose=verbose)
                return pd.DataFrame()
            
            self.log.write(" -Retrieved {} associations from API".format(len(associations_df)), verbose=verbose)
            
            # Debug: show available columns
            if verbose:
                self.log.write(" -Available columns: {}".format(list(associations_df.columns)[:10]), verbose=verbose)
            
            # Filter by p-value
            pval_col = None
            for col in ['p_value', 'pvalue', 'p']:
                if col in associations_df.columns:
                    pval_col = col
                    break
            
            if pval_col:
                associations_df[pval_col] = pd.to_numeric(associations_df[pval_col], errors='coerce')
                before_filter = len(associations_df)
                associations_df = associations_df[associations_df[pval_col] < sig_level].copy()
                self.log.write(" -Filtered to {} associations with p < {} (from {})".format(len(associations_df), sig_level, before_filter), verbose=verbose)
            else:
                self.log.write(" -Warning: No p-value column found, skipping p-value filter", verbose=verbose)
            
            if len(associations_df) == 0:
                return pd.DataFrame()
            
            # Extract CHR and POS from locations and build records
            records = []
            skipped_no_rsid = 0
            skipped_no_location = 0
            
            for idx, row in associations_df.iterrows():
                try:
                    # Extract variant rsID - try multiple column names
                    rsid = None
                    for col in ['variant_rsID', 'rsid', 'rs_id', 'snp_id']:
                        if col in row.index:
                            val = row[col]
                            # Handle Series/array values - extract scalar value
                            try:
                                if hasattr(val, '__len__') and not isinstance(val, str):
                                    if len(val) > 0:
                                        val = val.iloc[0] if hasattr(val, 'iloc') else val[0]
                                    else:
                                        continue
                                # Check if value is valid (not None, not NaN)
                                if val is not None and pd.notna(val):
                                    rsid = str(val)
                                    break
                            except (TypeError, ValueError, AttributeError, IndexError):
                                continue
                    
                    # Try extracting from snp_effect_allele if still no rsid
                    if not rsid and 'snp_effect_allele' in row.index:
                        val = row['snp_effect_allele']
                        # Handle Series/array values - extract scalar value
                        try:
                            if hasattr(val, '__len__') and not isinstance(val, str):
                                if len(val) > 0:
                                    val = val.iloc[0] if hasattr(val, 'iloc') else val[0]
                                else:
                                    val = None
                            if val is not None and pd.notna(val):
                                if isinstance(val, list) and len(val) > 0:
                                    allele_str = val[0]
                                    if isinstance(allele_str, str) and '-' in allele_str:
                                        rsid = allele_str.split('-')[0]
                                elif isinstance(val, str) and '-' in val:
                                    rsid = val.split('-')[0]
                        except (TypeError, ValueError, AttributeError, IndexError):
                            pass
                    
                    if not rsid:
                        skipped_no_rsid += 1
                        continue
                    
                    # Extract CHR and POS from locations
                    chr_val, pos_val = None, None
                    
                    # Try locations column
                    if 'locations' in row.index:
                        val = row['locations']
                        # Handle different data types that locations might be stored as
                        try:
                            # If it's already a string (e.g., "12:111803962"), use it directly
                            if isinstance(val, str):
                                chr_val, pos_val = self._extract_chr_pos_from_locations(val)
                            # If it's a list, get the first element
                            elif isinstance(val, list):
                                if len(val) > 0:
                                    # First element might be a string or another list
                                    first_elem = val[0]
                                    if isinstance(first_elem, str):
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(first_elem)
                                    elif isinstance(first_elem, list) and len(first_elem) > 0:
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(first_elem[0])
                                    else:
                                        # Try to extract from the list itself
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(val)
                                else:
                                    val = None
                            # If it's a pandas Series, get the value
                            elif hasattr(val, 'iloc'):
                                if len(val) > 0:
                                    series_val = val.iloc[0]
                                    if isinstance(series_val, (list, str)):
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(series_val)
                                    else:
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(val)
                                else:
                                    val = None
                            # If it's some other array-like object
                            elif hasattr(val, '__len__') and not isinstance(val, str):
                                if len(val) > 0:
                                    array_val = val[0] if hasattr(val, '__getitem__') else None
                                    if array_val is not None:
                                        chr_val, pos_val = self._extract_chr_pos_from_locations(array_val)
                                else:
                                    val = None
                            
                            # If we still don't have chr/pos, try the original value
                            if (not chr_val or not pos_val) and val is not None and pd.notna(val):
                                chr_val, pos_val = self._extract_chr_pos_from_locations(val)
                        except (TypeError, ValueError, AttributeError, IndexError) as e:
                            pass
                    
                    if not chr_val or not pos_val:
                        skipped_no_location += 1
                        continue
                    
                    # Extract other metadata
                    pval = None
                    for col in ['p_value', 'pvalue', 'p']:
                        if col in row.index:
                            val = row[col]
                            try:
                                # Handle Series/array values
                                if hasattr(val, '__len__') and not isinstance(val, str):
                                    if len(val) > 0:
                                        val = val.iloc[0] if hasattr(val, 'iloc') else val[0]
                                    else:
                                        continue
                                if val is not None and pd.notna(val):
                                    pval = val
                                    break
                            except (TypeError, ValueError, AttributeError, IndexError):
                                continue
                    
                    beta = row.get('beta', None) or row.get('betaNum', None)
                    se = row.get('standard_error', None) or row.get('se', None) or row.get('standardError', None)
                    or_val = row.get('or_per_copy_num', None) or row.get('or', None) or row.get('orPerCopyNum', None)
                    
                    # Extract trait information
                    trait = None
                    for col in ['reported_trait', 'trait', 'disease_trait']:
                        if col in row.index:
                            val = row[col]
                            # Handle Series/array values - extract scalar value
                            try:
                                if hasattr(val, '__len__') and not isinstance(val, str):
                                    if len(val) > 0:
                                        val = val.iloc[0] if hasattr(val, 'iloc') else val[0]
                                    else:
                                        continue
                                if val is not None and pd.notna(val):
                                    if isinstance(val, list):
                                        trait = ', '.join([str(t) for t in val])
                                    else:
                                        trait = str(val)
                                    break
                            except (TypeError, ValueError, AttributeError, IndexError):
                                continue
                    
                    # Extract gene information
                    gene = None
                    for col in ['mapped_genes', 'gene', 'genes', 'mappedGenes']:
                        if col in row.index:
                            val = row[col]
                            # Handle Series/array values - extract scalar value
                            try:
                                if hasattr(val, '__len__') and not isinstance(val, str):
                                    if len(val) > 0:
                                        val = val.iloc[0] if hasattr(val, 'iloc') else val[0]
                                    else:
                                        continue
                                if val is not None and pd.notna(val):
                                    if isinstance(val, list):
                                        gene = ', '.join([str(g) for g in val])
                                    else:
                                        gene = str(val)
                                    break
                            except (TypeError, ValueError, AttributeError, IndexError):
                                continue
                    
                    # Extract study information
                    pubmed_id = row.get('pubmed_id', None) or row.get('pubmedId', None)
                    author = row.get('first_author', None) or row.get('author', None) or row.get('firstAuthor', None)
                    study_title = row.get('study_accession', None) or row.get('accession_id', None) or row.get('studyAccession', None)
                
                    records.append({
                        'SNPID': rsid,
                        'CHR': chr_val,
                        'POS': pos_val,
                        'P': pval,
                        'BETA': beta,
                        'SE': se,
                        'OR': or_val,
                        'TRAIT': trait,
                        'REPORT_GENENAME': gene,
                        'PUBMEDID': pubmed_id,
                        'AUTHOR': author,
                        'STUDY': study_title
                    })
                except Exception as e:
                    self.log.write(" -Warning: Error processing association row {}: {}".format(idx, str(e)), verbose=verbose)
                    continue
            
            if skipped_no_rsid > 0:
                self.log.write(" -Skipped {} associations without valid rsID".format(skipped_no_rsid), verbose=verbose)
            if skipped_no_location > 0:
                self.log.write(" -Skipped {} associations without valid CHR/POS".format(skipped_no_location), verbose=verbose)
            
            if len(records) == 0:
                return pd.DataFrame()
            
            result_df = pd.DataFrame(records)
            self.log.write(" -Successfully extracted {} variants with CHR/POS".format(len(result_df)), verbose=verbose)
            return result_df
            
        except Exception as e:
            self.log.write(" -Error retrieving data from GWAS Catalog API v2: {}. Falling back to old API...".format(str(e)), verbose=verbose)
            # Fallback to old API
            try:
                from gwaslab.util.util_ex_gwascatalog import gwascatalog_trait
                known_Sumstats = gwascatalog_trait(efo, source="NCBI", sig_level=sig_level, 
                                                   use_cache=True, cache_dir="./", 
                                                   verbose=verbose, log=self.log)
                if known_Sumstats and hasattr(known_Sumstats, 'data') and len(known_Sumstats.data) > 0:
                    self.log.write(" -Retrieved {} variants using old API v1".format(len(known_Sumstats.data)), verbose=verbose)
                    return known_Sumstats.data.copy()
            except Exception as e2:
                self.log.write(" -Fallback to old API also failed: {}".format(str(e2)), verbose=verbose)
            
            return pd.DataFrame()


# ============================================================================
# Convenience Functions
# ============================================================================

def get_studies(**kwargs) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
    """Convenience function to get studies."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_studies(**kwargs)


def get_associations(**kwargs) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
    """Convenience function to get associations."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_associations(**kwargs)


def get_variants(**kwargs) -> Union[pd.DataFrame, Dict[str, Any], List[Dict[str, Any]]]:
    """Convenience function to get variants."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_variants(**kwargs)


def search_studies(query: str, **kwargs) -> pd.DataFrame:
    """Convenience function to search studies by trait."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_studies(efo_trait=query, **kwargs)


def search_associations(query: str, **kwargs) -> pd.DataFrame:
    """Convenience function to search associations by trait."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_associations(efo_trait=query, **kwargs)


def search_variants(query: str, **kwargs) -> pd.DataFrame:
    """Convenience function to search variants by rsID."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.get_variants(rs_id=query, **kwargs)


def search_traits(query: str, **kwargs) -> pd.DataFrame:
    """Convenience function to search traits by free text."""
    verbose = kwargs.pop("verbose", True)
    client = GWASCatalogClient(verbose=verbose)
    return client.search_traits(query, **kwargs)
