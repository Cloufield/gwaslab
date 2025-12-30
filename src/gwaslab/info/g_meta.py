from typing import TYPE_CHECKING, Dict, Any, Optional
from gwaslab.info.g_version import gwaslab_info
from datetime import datetime
import numpy as np
from gwaslab.util.util_in_filter_value import _infer_build
from gwaslab.info.g_Log import Log
import time
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _init_meta(object: str = "Sumstats") -> Dict[str, Any]:
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

def _append_meta_record(old: str, new: str) -> str:
    if old == "Unknown" or old== "Unchecked":
        return new
    else:
        return "{}, {}".format(old, new)

def _update_step_status(status_section: Optional[Dict[str, Any]], step_name: str, now: str, performed: bool, params: Dict[str, Any]) -> None:
    if status_section is None:
        return
    step = status_section.get(step_name)
    if step is None:
        return
    step["performed"] = bool(performed)
    if performed:
        step["last_executed"] = now
    step["parameters_used"] = params if params is not None else {}

def _update_qc_step(self: 'Sumstats', step_name: str, params: Dict[str, Any], performed: bool = True) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    _update_step_status(qc_status, step_name, now, performed, params)

def _update_harmonize_step(self: 'Sumstats', step_name: str, params: Dict[str, Any], performed: bool = True) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    harm_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("harmonize", {})
    _update_step_status(harm_status, step_name, now, performed, params)

def _check_sumstats_qc_status(self: 'Sumstats') -> Dict[str, Any]:
    """
    Check the QC and harmonization status of the sumstats.
    
    Returns
    -------
    dict
        Dictionary containing QC and harmonization status information with keys:
        - "basic_check": Basic QC check status
        - "harmonize": Harmonization status
        - "qc_and_harmonization_status": Overall QC and harmonization status
    """
    return {
        "basic_check": self.meta["gwaslab"].get("basic_check", {}),
        "harmonize": self.meta["gwaslab"].get("harmonize", {}),
        "qc_and_harmonization_status": self.meta["gwaslab"].get("qc_and_harmonization_status", {})
    }

def _set_qc_status(self: 'Sumstats', args: Dict[str, Any]) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    self.meta["gwaslab"]["basic_check"]["performed"] = True
    self.meta["gwaslab"]["basic_check"]["last_executed"] = now
    args_to_save = {k: v for k, v in args.items() if k != "self"}
    self.meta["gwaslab"]["basic_check"]["parameters_used"] = args_to_save

    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    if qc_status:
        _update_step_status(qc_status, "id", now, True, args_to_save.get("fix_id_kwargs", {}))
        _update_step_status(qc_status, "chr", now, True, args_to_save.get("fix_chr_kwargs", {}))
        _update_step_status(qc_status, "pos", now, True, args_to_save.get("fix_pos_kwargs", {}))
        _update_step_status(qc_status, "allele", now, True, args_to_save.get("fix_allele_kwargs", {}))
        _update_step_status(qc_status, "sanity", now, True, args_to_save.get("sanity_check_stats_kwargs", {}))
        _update_step_status(qc_status, "consistency", now, True, args_to_save.get("consistency_check_kwargs", {}))
        _update_step_status(qc_status, "normalize", now, bool(args_to_save.get("normalize", False)), args_to_save.get("normalize_allele_kwargs", {}))
        _update_step_status(qc_status, "remove_dup", now, bool(args_to_save.get("remove_dup", False)), args_to_save.get("remove_dup_kwargs", {}))
        _update_step_status(qc_status, "sort_coord", now, True, {})
        _update_step_status(qc_status, "sort_column", now, True, {})

def _set_harmonization_status(self: 'Sumstats', args: Dict[str, Any]) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    self.meta["gwaslab"]["harmonize"]["performed"]= True
    self.meta["gwaslab"]["harmonize"]["last_executed"] = now
    args_to_save = {k: v for k, v in args.items() if k != "self"}
    self.meta["gwaslab"]["harmonize"]["parameters_used"] = args_to_save

    harm_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("harmonize", {})
    if harm_status:
        if args_to_save.get("ref_seq", None) is not None:
            _update_step_status(harm_status, "check_ref", now, True, args_to_save.get("check_ref_kwargs", {}))
            _update_step_status(harm_status, "flip_allele_stats", now, True, args_to_save.get("flip_allele_stats_kwargs", {}))

        if args_to_save.get("ref_infer", None) is not None:
            _update_step_status(harm_status, "infer_strand", now, True, args_to_save.get("infer_strand_kwargs", {}))

        if args_to_save.get("ref_rsid_tsv", None) is not None or args_to_save.get("ref_rsid_vcf", None) is not None:
            _update_step_status(harm_status, "assign_rsid", now, True, args_to_save.get("assign_rsid_kwargs", {}))

        _update_step_status(harm_status, "liftover", now, False, args_to_save.get("liftover_kwargs", {}))

def _update_meta(meta: Dict[str, Any], sumstats: pd.DataFrame, object: str = "Sumstats", log: Log = Log(), verbose: bool = True) -> Dict[str, Any]:  
    """
    Update Sumstats Object meta info based on the statistics of the current sumstats.
    Including information on variants, samples.
    """
    meta["gwaslab"]["variants"]["variant_number"] = len(sumstats)
    
    if "CHR" in sumstats.columns:
        meta["gwaslab"]["variants"]["number_of_chromosomes"] = len(sumstats["CHR"].unique())
    
    if  meta["gwaslab"]["gwaslab_object"]=="gwaslab.Sumstats":      
        if "P" in sumstats.columns:
            meta["gwaslab"]["variants"]["min_P"]=np.nanmin(sumstats["P"])
        if "EAF" in sumstats.columns:
            meta["gwaslab"]["variants"]["min_minor_allele_freq"]=min (np.min(sumstats["EAF"]) , 1- np.max(sumstats["EAF"]))
        if "N" in sumstats.columns:
            meta["gwaslab"]["samples"]["sample_size"] = int(sumstats["N"].max())
            meta["gwaslab"]["samples"]["sample_size_median"] = sumstats["N"].median()
            meta["gwaslab"]["samples"]["sample_size_min"] = int(sumstats["N"].min())
    
    if  meta["gwaslab"]["gwaslab_object"]=="gwaslab.SumstatsMulti" or meta["gwaslab"]["gwaslab_object"]=="gwaslab.SumstatsPair":   
        nstudy = meta["gwaslab"]['number_of_studies']
        for i in range(nstudy):
            i_form_1 = i + 1
            meta["gwaslab"]["variants"][i_form_1]=dict()
            meta["gwaslab"]["samples"][i_form_1] =dict()

            if "P_{}".format(i_form_1) in sumstats.columns:
                p = "P_{}".format(i_form_1)
                
                meta["gwaslab"]["variants"][i_form_1]["min_P"]= np.nanmin(sumstats[p])
            if "N_{}".format(i_form_1) in sumstats.columns:
                n = "N_{}".format(i_form_1)
                meta["gwaslab"]["samples"][i_form_1]["sample_size"] = int(sumstats[n].max())
                meta["gwaslab"]["samples"][i_form_1]["sample_size_median"] = sumstats[n].median()
                meta["gwaslab"]["samples"][i_form_1]["sample_size_min"] = int(sumstats[n].min())
            if "EAF_{}".format(i_form_1) in sumstats.columns:
                eaf="EAF_{}".format(i_form_1)
                meta["gwaslab"]["variants"][i_form_1]["min_minor_allele_freq"]=min (np.min(sumstats[eaf]) , 1- np.max(sumstats[eaf]))

    if meta["gwaslab"]["genome_build"] == "99":
        _, meta["gwaslab"]["genome_build"] = _infer_build(sumstats, change_status=False, log=log, verbose=verbose)
    
    meta["date_last_modified"] = str(time.strftime('%Y/%m/%d'))
    
    return meta
