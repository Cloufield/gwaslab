# Sumstats — Downstream

Lead/novel loci, associations, LDSC, clumping, finemapping, and PRS.

!!! note "Extension runners not in API Reference"
    `run_*` methods (SuSiE, PRS-CS, MAGMA, scDRS) are omitted from this page until ready for publication. They remain available on `Sumstats` objects.

::: gwaslab.Sumstats
    options:
      show_root_heading: false
      show_root_full_path: false
      heading_level: 2
      inherited_members: false
      members:
        - get_lead
        - get_top
        - get_novel
        - get_density
        - get_associations
        - check_cis
        - check_novel_set
        - check_cs_overlap
        - anno_gene
        - get_per_snp_r2
        - get_ess
        - get_gc
        - infer_ancestry
        - abf_finemapping
        - get_cs_lead
        - read_pipcs
        - clump
        - calculate_prs
        - estimate_h2_by_ldsc
        - estimate_rg_by_ldsc
        - estimate_h2_cts_by_ldsc
        - estimate_partitioned_h2_by_ldsc
        - calculate_ld_matrix
        - extract_ld_matrix
        - get_ld_matrix_from_vcf
      filters:
        - "!^_[^_]"
