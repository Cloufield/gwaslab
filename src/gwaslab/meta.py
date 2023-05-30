from gwaslab.version import gwaslab_info

def init_meta():
    metadata = {"gwaslab":{
                        "gwaslab_version": gwaslab_info()["version"],
                        "study_name":"Sumstats_1",
                        "study_type":"Unknown",
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
                            "ref_infer":"Unknown"
                        }
                    },
                     "genotyping_technology":"Unknown", 
                     "gwas_id":"Unknown", 
                     "samples":{
                          "sample_size":"Unknown", 
                          "sample_ancestry":"European", 
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
    return metadata.copy()