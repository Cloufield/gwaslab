from typing import TYPE_CHECKING, Dict, Any, Optional
from gwaslab.info.g_version import gwaslab_info
from datetime import datetime
import numpy as np
from gwaslab.util.util_in_filter_value import _match_build_from_hapmap3
from gwaslab.qc.qc_build import _build_to_assembly
from gwaslab.info.g_Log import Log
import time
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _init_meta(object: str = "Sumstats") -> Dict[str, Any]:
    """Initialize the meta information of the Sumstats Object.
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
                    "genome_build_confidence":"unknown",
                    "genome_build_source":"Unknown",
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
                            "rsid_to_chrpos": {"performed": False, "last_executed": "", "parameters_used": {}},
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

def _update_step_status(status_section: Optional[Dict[str, Any]], step_name: str, now: str, performed: bool, params: Dict[str, Any], log: Optional[Log] = None) -> None:
    if status_section is None:
        return
    step = status_section.get(step_name)
    if step is None:
        if log is not None:
            log.warning(f"Unknown metadata step '{step_name}' — schema may be out of date.", verbose=False)
        return
    step["performed"] = bool(performed)
    if performed:
        step["last_executed"] = now
    step["parameters_used"] = params if params is not None else {}

def _update_qc_step(self: 'Sumstats', step_name: str, params: Dict[str, Any], performed: bool = True) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    _update_step_status(qc_status, step_name, now, performed, params, log=getattr(self, "log", None))

def _update_harmonize_step(self: 'Sumstats', step_name: str, params: Dict[str, Any], performed: bool = True) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    harm_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("harmonize", {})
    _update_step_status(harm_status, step_name, now, performed, params, log=getattr(self, "log", None))

def _check_sumstats_qc_status(self: 'Sumstats') -> Dict[str, Any]:
    """Check the QC and harmonization status of the sumstats.
    
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

def _set_qc_status(self: 'Sumstats', args: Dict[str, Any], qc_report: Optional[Any] = None) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    self.meta["gwaslab"]["basic_check"]["performed"] = True
    self.meta["gwaslab"]["basic_check"]["last_executed"] = now
    args_to_save = {k: v for k, v in args.items() if k != "self"}
    self.meta["gwaslab"]["basic_check"]["parameters_used"] = args_to_save

    qc_status = self.meta["gwaslab"].get("qc_and_harmonization_status", {}).get("qc", {})
    if not qc_status:
        return

    step_param_map = {
        "id": "fix_id_kwargs",
        "chr": "fix_chr_kwargs",
        "pos": "fix_pos_kwargs",
        "allele": "fix_allele_kwargs",
        "sanity": "sanity_check_stats_kwargs",
        "consistency": "consistency_check_kwargs",
        "normalize": "normalize_allele_kwargs",
        "remove_dup": "remove_dup_kwargs",
        "sort_coord": None,
        "sort_column": None,
    }

    if qc_report is not None and hasattr(qc_report, "steps"):
        for step_name, ran in qc_report.steps.items():
            param_key = step_param_map.get(step_name)
            params = args_to_save.get(param_key, {}) if param_key else {}
            if step_name == "normalize":
                params = args_to_save.get("normalize_allele_kwargs", {})
            elif step_name == "remove_dup":
                params = args_to_save.get("remove_dup_kwargs", {})
            _update_step_status(qc_status, step_name, now, ran, params if ran else {})
    else:
        for step_name, param_key in step_param_map.items():
            if step_name in ("normalize", "remove_dup"):
                ran = bool(args_to_save.get(step_name if step_name != "normalize" else "normalize", False))
                if step_name == "normalize":
                    ran = bool(args_to_save.get("normalize", False))
                elif step_name == "remove_dup":
                    ran = bool(args_to_save.get("remove_dup", False))
            elif step_name == "consistency":
                ran = bool(args_to_save.get("include_consistency", True))
            else:
                ran = True
            params = args_to_save.get(param_key, {}) if param_key else {}
            _update_step_status(qc_status, step_name, now, ran, params if ran else {})

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
    """Update Sumstats Object meta info based on the statistics of the current sumstats.
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

    if meta["gwaslab"]["genome_build"] == "99" and "CHR" in sumstats.columns and "POS" in sumstats.columns:
        inferred, confidence, _, _ = _match_build_from_hapmap3(
            sumstats, log=log, verbose=verbose
        )
        if inferred not in ("Unknown", "99"):
            species = meta["gwaslab"].get("species", "homo sapiens")
            meta["gwaslab"]["genome_build"] = inferred
            meta["genome_assembly"] = _build_to_assembly(species, inferred)
            meta["gwaslab"]["genome_build_confidence"] = confidence
            meta["gwaslab"]["genome_build_source"] = "update_meta"
    
    meta["date_last_modified"] = str(time.strftime('%Y/%m/%d'))
    
    return meta


def _validate_meta(self: 'Sumstats') -> Dict[str, Any]:
    """Validate internal metadata consistency (build fields, QC flags).

Returns
-------
dict
    Structured report with keys ``valid`` (bool), ``issues`` (list), ``build`` (dict).
"""
    from gwaslab.qc.qc_build import _build_to_assembly

    issues: list = []
    gwaslab_meta = self.meta.get("gwaslab", {})
    build_code = getattr(self, "_build", gwaslab_meta.get("genome_build", "99"))
    meta_build = gwaslab_meta.get("genome_build", "99")
    assembly = self.meta.get("genome_assembly", "Unknown")
    species = gwaslab_meta.get("species", "homo sapiens")
    expected_assembly = _build_to_assembly(species, meta_build)

    if str(build_code) != str(meta_build):
        issues.append(
            f"Build mismatch: _build={build_code!r} vs meta genome_build={meta_build!r}"
        )
    if assembly != expected_assembly and not (assembly == "Unknown" and meta_build == "99"):
        issues.append(
            f"genome_assembly={assembly!r} does not match expected {expected_assembly!r} for build {meta_build}"
        )

    basic = gwaslab_meta.get("basic_check", {})
    qc_steps = gwaslab_meta.get("qc_and_harmonization_status", {}).get("qc", {})
    if basic.get("performed") and qc_steps:
        any_qc = any(step.get("performed") for step in qc_steps.values())
        if not any_qc:
            issues.append("basic_check.performed=True but no QC sub-steps are marked performed")

    confidence = gwaslab_meta.get("genome_build_confidence", "unknown")
    source = gwaslab_meta.get("genome_build_source", "Unknown")

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "build": {
            "code": meta_build,
            "assembly": assembly,
            "confidence": confidence,
            "source": source,
        },
    }
