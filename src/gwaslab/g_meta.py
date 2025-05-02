from gwaslab.g_version import gwaslab_info

def _init_meta(object="Sumstats"):
    metadata_ssf ={
                        "genotyping_technology":"Unknown", 
                        "gwas_id":"Unknown", 
                        "samples":{
                            "sample_size":"Unknown", 
                            "sample_ancestry":"Unknown", 
                            "ancestry_method":"self-reported|genetically determined", 
                        } ,
                        "trait_description":"Unknown", 
                        "minor_allele_freq_lower_limit":"Unknown", 
                        "data_file_name":"Unknown", 
                        "file_type":"Unknown", 
                        "data_file_md5sum":"Unknown", 
                        "is_harmonised":"Unchecked", 
                        "is_sorted":"Unchecked", 
                        "genome_assembly":"Unknown",
                        "date_last_modified":"Unknown", 
                        "coordinate_system":"1-based",
                        "sex": "M|F|combined"
                    }
    metadata_multi ={
                        "genome_assembly":"Unknown",
                        "date_last_modified":"Unknown", 
                        "coordinate_system":"1-based"
                    }
   
    # Sumstats
    if object=="Sumstats":
        metadata = {"gwaslab":{
                    "gwaslab_version": gwaslab_info()["version"],
                    "gwaslab_object":"gwaslab.Sumstats",
                    "study_name":"Sumstats1",
                    "study_type":"Unknown",
                    "species":"homo sapiens",
                    "genome_build":"99",
                    "sample_prevalence":"Unknown",
                    "inferred_ancestry":"Unknown",
                    "population_prevalence":"Unknown",
                    "variants":{
                        "variant_number":"Unknown",
                        "min_P":"Unknown",
                        "number_of_chromosomes":"Unknown",
                    },
                    "samples":{
                        "sample_size":"Unknown",
                        "sample_size_case":"Unknown",
                        "sample_size_control":"Unknown",
                        "sample_size_median":"Unknown",
                        "sample_size_min":"Unknown",
                    },
                    "references":{
                        "ref_rsid_tsv":"Unknown",
                        "ref_rsid_vcf":"Unknown",
                        "ref_seq":"Unknown",
                        "ref_infer":"Unknown",
                        "ref_infer_af":"Unknown",
                        "ref_infer_daf":"Unknown",
                        "ref_rsid_to_chrpos_tsv":"Unknown",
                        "ref_rsid_to_chrpos_vcf":"Unknown"
                    }}}
        metadata |= metadata_ssf
    
    # SumstatsPair
    elif object=="SumstatsPair":
        metadata = {"gwaslab":{
                                    "gwaslab_version": gwaslab_info()["version"],
                                    "gwaslab_object":"gwaslab.SumstatsPair",
                                    "group_name":"Group1",
                                    "species":"homo sapiens",
                                    "genome_build":"99",
                                    "variants":{
                                        "variant_number":"Unknown",
                                        "min_P":"Unknown",
                                        "number_of_chromosomes":"Unknown",
                                    },
                                    "samples":{
                                        "sample_size":"Unknown",
                                        "sample_size_case":"Unknown",
                                        "sample_size_control":"Unknown",
                                        "sample_size_median":"Unknown",
                                        "sample_size_min":"Unknown",
                                    },
                                    "references":{
                                        "ref_rsid_tsv":"Unknown",
                                        "ref_rsid_vcf":"Unknown",
                                        "ref_seq":"Unknown",
                                        "ref_infer":"Unknown",
                                        "ref_infer_af":"Unknown",
                                        "ref_infer_daf":"Unknown",
                                        "ref_rsid_to_chrpos_tsv":"Unknown",
                                        "ref_rsid_to_chrpos_vcf":"Unknown"
                                    }}}
        metadata |= metadata_multi
    
    # SumstatsMulti
    elif object=="SumstatsMulti":
        metadata = {"gwaslab":{
                            "gwaslab_version": gwaslab_info()["version"],
                            "gwaslab_object":"gwaslab.SumstatsMulti",
                            "group_name":"Group1",
                            "species":"homo sapiens",
                            "genome_build":"99",
                            "variants":{
                                "variant_number":"Unknown",
                                "min_P":"Unknown",
                                "number_of_chromosomes":"Unknown",
                            },
                            "samples":{
                                "sample_size":"Unknown",
                                "sample_size_case":"Unknown",
                                "sample_size_control":"Unknown",
                                "sample_size_median":"Unknown",
                                "sample_size_min":"Unknown",
                            },
                            "references":{
                                "ref_rsid_tsv":"Unknown",
                                "ref_rsid_vcf":"Unknown",
                                "ref_seq":"Unknown",
                                "ref_infer":"Unknown",
                                "ref_infer_af":"Unknown",
                                "ref_infer_daf":"Unknown",
                                "ref_rsid_to_chrpos_tsv":"Unknown",
                                "ref_rsid_to_chrpos_vcf":"Unknown"
                            }}}
        metadata |= metadata_multi
    return metadata.copy()

def _append_meta_record(old, new):
    if old == "Unknown" or old== "Unchecked":
        return new
    else:
        return "{}, {}".format(old, new)
    