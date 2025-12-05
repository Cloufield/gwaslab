from gwaslab.g_version import gwaslab_info
from datetime import datetime


def _init_meta(object="Sumstats"):
    """
    Initialize the meta information of the Sumstats Object. 
    """
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
                    },
                    "basic_check": {
                            "performed": False,
                            "last_executed": "",
                            "parameters_used": {}
                    },
                    "harmonize": {
                        "performed": False,
                        "last_executed": "",
                        "parameters_used": {}
                    },
                    "qc_and_harmonization_status": {
                        "qc": {
                            "id": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "chr": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "pos": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "allele": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "sanity": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "consistency": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "normalize": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "remove_dup": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "sort_coord": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "sort_column": {"performed": False, "last_executed": "", "parameters_used": {}}
                        },
                        "harmonize": {
                            "check_ref": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "flip_allele_stats": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "infer_strand": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "assign_rsid": {"performed": False, "last_executed": "", "parameters_used": {}},
                            "liftover": {"performed": False, "last_executed": "", "parameters_used": {}}
                        }
                    }
                    }}
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
                                    }                             
                                    }}
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

def _update_step_status(status_section, step_name, now, performed, params):
    if status_section is None:
        return
    step = status_section.get(step_name)
    if step is None:
        return
    step["performed"] = bool(performed)
    if performed:
        step["last_executed"] = now
    step["parameters_used"] = params if params is not None else {}

def _update_qc_step(self, step_name, params, performed=True):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    _update_step_status(qc_status, step_name, now, performed, params)

def _update_harmonize_step(self, step_name, params, performed=True):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    harm_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("harmonize", {})
    _update_step_status(harm_status, step_name, now, performed, params)

def _check_sumstats_qc_status(self):
    return {
        "basic_check": self.meta["gwaslab"].get("basic_check", {}),
        "harmonize": self.meta["gwaslab"].get("harmonize", {}),
        "qc_and_harmonization_status": self.meta["gwaslab"].get("qc_and_harmonization_status", {})
    }

def _set_qc_status(self, args):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    self.meta["gwaslab"]["basic_check"]["performed"] = True
    self.meta["gwaslab"]["basic_check"]["last_executed"] = now
    args_to_save = {k: v for k, v in args.items() if k != "self"}
    self.meta["gwaslab"]["basic_check"]["parameters_used"] = args_to_save

    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    if qc_status:
        _update_step_status(qc_status, "id", now, True, args_to_save.get("fixid_args", {}))
        _update_step_status(qc_status, "chr", now, True, args_to_save.get("fixchr_args", {}))
        _update_step_status(qc_status, "pos", now, True, args_to_save.get("fixpos_args", {}))
        _update_step_status(qc_status, "allele", now, True, args_to_save.get("fixallele_args", {}))
        _update_step_status(qc_status, "sanity", now, True, args_to_save.get("sanitycheckstats_args", {}))
        _update_step_status(qc_status, "consistency", now, True, args_to_save.get("consistencycheck_args", {}))
        _update_step_status(qc_status, "normalize", now, bool(args_to_save.get("normalize", False)), args_to_save.get("normalizeallele_args", {}))
        _update_step_status(qc_status, "remove_dup", now, bool(args_to_save.get("remove_dup", False)), args_to_save.get("removedup_args", {}))
        _update_step_status(qc_status, "sort_coord", now, True, {})
        _update_step_status(qc_status, "sort_column", now, True, {})

def _set_harmonization_status(self, args):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    self.meta["gwaslab"]["harmonize"]["performed"]= True
    self.meta["gwaslab"]["harmonize"]["last_executed"] = now
    args_to_save = {k: v for k, v in args.items() if k != "self"}
    self.meta["gwaslab"]["harmonize"]["parameters_used"] = args_to_save

    harm_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("harmonize", {})
    if harm_status:
        if args_to_save.get("ref_seq", None) is not None:
            _update_step_status(harm_status, "check_ref", now, True, args_to_save.get("checkref_args", {}))
            _update_step_status(harm_status, "flip_allele_stats", now, True, args_to_save.get("flipallelestats_args", {}))

        if args_to_save.get("ref_infer", None) is not None:
            _update_step_status(harm_status, "infer_strand", now, True, args_to_save.get("inferstrand_args", {}))

        if args_to_save.get("ref_rsid_tsv", None) is not None or args_to_save.get("ref_rsid_vcf", None) is not None:
            _update_step_status(harm_status, "assign_rsid", now, True, args_to_save.get("assignrsid_args", {}))

        _update_step_status(harm_status, "liftover", now, False, args_to_save.get("liftover_args", {}))
