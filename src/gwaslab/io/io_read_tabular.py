import gzip
from typing import Any, Dict, Mapping, Optional, Union

import pandas as pd
from gwaslab.bd.bd_common_data import get_formats_list
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_common_data import get_format_dict


def _pre_rename_dtype_map(
    meta_data: Mapping[str, Any], dtypes: Mapping[str, Any]
) -> Dict[Union[str, int], Any]:
    """Map formatbook ``format_datatype`` keys to labels pandas used before rename."""
    if "format_header" in meta_data:
        fh = meta_data["format_header"]
        if fh is None or fh is False:
            out: Dict[Union[str, int], Any] = {}
            for k, v in dtypes.items():
                try:
                    out[int(k)] = v
                except (TypeError, ValueError):
                    out[str(k)] = v
            return out
    return dict(dtypes)


def format_comment_leading_skip_count(meta_data: Optional[Mapping[str, Any]], path: str) -> int:
    """
    Number of leading file lines to skip for formatbook ``format_comment``.

    ``format_comment='#'`` skips only consecutive ``##`` meta lines (VCF header /
    PLINK comments), preserving the ``#CHROM`` column header row.
    """
    if not meta_data:
        return 0
    fc = meta_data.get("format_comment")
    if not isinstance(fc, str):
        return 0
    if fc == "#":
        return _count_leading_lines_with_prefix(path, "##")
    if len(fc) > 1:
        return _count_leading_lines_with_prefix(path, fc)
    return 0


def _merge_skiprows_into_readargs(readargs: Dict[str, Any], new_rows: list) -> None:
    existing = readargs.get("skiprows")
    skip = set(new_rows)
    if existing is not None:
        if isinstance(existing, int):
            skip.update(range(existing))
        else:
            skip.update(existing)
    readargs["skiprows"] = sorted(skip)


def apply_format_comment_readargs(
    meta_data: Optional[Mapping[str, Any]],
    readargs: Dict[str, Any],
    path: Optional[str] = None,
    user_kwargs: Optional[Mapping[str, Any]] = None,
) -> None:
    """
    Apply formatbook ``format_comment`` to pandas read kwargs.

    ``format_comment='#'`` must not become ``comment='#'`` (that drops ``#CHROM``).
    """
    user_kwargs = user_kwargs or {}
    if not meta_data:
        return
    fc = meta_data.get("format_comment")
    if fc is None or not isinstance(fc, str):
        return

    if fc == "#":
        readargs.pop("comment", None)
        return

    if len(fc) > 1 and path and "skiprows" not in user_kwargs:
        skip_n = _count_leading_lines_with_prefix(path, fc)
        if skip_n:
            _merge_skiprows_into_readargs(readargs, list(range(skip_n)))
    elif len(fc) == 1 and "comment" not in user_kwargs:
        readargs["comment"] = fc


def build_path_skiprows(inpath: str, meta_data: Optional[Mapping[str, Any]]) -> Optional[list]:
    """Leading skip row indices from ``format_comment`` (excludes vcf.gz handling)."""
    n = format_comment_leading_skip_count(meta_data, inpath) if inpath else 0
    return list(range(n)) if n else None


def _count_leading_lines_with_prefix(path: str, prefix: str) -> int:
    """Count initial lines starting with ``prefix`` (e.g. skip VCF/PLINK2 ``##`` meta)."""
    n = 0
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(prefix):
                    n += 1
                else:
                    break
    else:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(prefix):
                    n += 1
                else:
                    break
    return n


def _multiline_header_skiprows(
    meta_data: Mapping[str, Any],
    load_kwargs_dict: Dict[str, Any],
    user_kwargs: Mapping[str, Any],
) -> None:
    """
    Skip trailing header rows after the column-name row (formatbook ``format_header_lines``).

    Line 0 of the file (after any leading ``skiprows``) is the column header; lines
    ``lead + 1`` … ``lead + format_header_lines - 1`` are additional header lines to skip.
    """
    if "skiprows" in user_kwargs:
        return
    raw_nl = meta_data.get("format_header_lines", 1)
    try:
        n_header_lines = int(raw_nl)
    except (TypeError, ValueError):
        return
    if n_header_lines <= 1:
        return
    if "format_header" not in meta_data:
        return
    fh = meta_data["format_header"]
    if fh is None or fh is False:
        return

    def _header_row_index(skiprows_val: Any) -> Optional[int]:
        if skiprows_val is None:
            return 0
        if isinstance(skiprows_val, int):
            return skiprows_val
        if isinstance(skiprows_val, list):
            if not skiprows_val:
                return 0
            if skiprows_val == list(range(len(skiprows_val))):
                return len(skiprows_val)
        return None

    existing = load_kwargs_dict.get("skiprows")
    lead = _header_row_index(existing)
    if lead is None:
        return
    extra = list(range(lead + 1, lead + n_header_lines))
    if not extra:
        return
    if existing is None:
        load_kwargs_dict["skiprows"] = extra
    elif isinstance(existing, int):
        load_kwargs_dict["skiprows"] = list(range(existing)) + extra
    else:
        load_kwargs_dict["skiprows"] = list(existing) + extra


def _read_tabular(path: str, fmt: str, **kwargs: Any) -> pd.DataFrame:
    
    # default
    load_kwargs_dict = {"sep":"\t",
                      "header":None}
    
    # if specified by user
    if len(kwargs)>0:
        load_kwargs_dict = kwargs
    
    # load format
    meta_data, rename_dictionary = get_format_dict(fmt)
    
    if "format_separator" in meta_data and "sep" not in kwargs:
        load_kwargs_dict["sep"] = meta_data["format_separator"]
    
    apply_format_comment_readargs(meta_data, load_kwargs_dict, path, kwargs)
    fc_skip = build_path_skiprows(path, meta_data)
    if fc_skip and "skiprows" not in kwargs:
        _merge_skiprows_into_readargs(load_kwargs_dict, fc_skip)

    if "format_header" in meta_data and "header" not in kwargs:
        if meta_data["format_header"] is True:
            load_kwargs_dict["header"] = "infer"
        elif meta_data["format_header"] is False:
            load_kwargs_dict["header"] = None
        else:
            load_kwargs_dict["header"] = meta_data["format_header"]

    if "format_na" in meta_data and "na_values" not in kwargs:
        if  meta_data["format_na"] is not None:    
            load_kwargs_dict["na_values"] = meta_data["format_na"]

    _multiline_header_skiprows(meta_data, load_kwargs_dict, kwargs)

    #######################################################################################
    df = pd.read_csv(path, **load_kwargs_dict)
    #######################################################################################

    # format_datatype keys = on-disk names before rename; only cast columns that exist
    if "format_datatype" in meta_data:
        dtype_map = _pre_rename_dtype_map(meta_data, meta_data["format_datatype"])
        dtype_map = {k: v for k, v in dtype_map.items() if k in df.columns}
        if dtype_map:
            df = df.astype(dtype_map)

    # rename columns
    # False or None => no header row was read; pandas columns are 0,1,2,... and format_dict
    # keys are "0","1",... — convert to int so rename matches. True => file supplies names;
    # rename_dictionary keys match those header strings.
    if "format_header" in meta_data:
        fh = meta_data["format_header"]
        if fh is None or fh is False:
            num_to_name = {int(k): v for k, v in rename_dictionary.items()}
            df = df.rename(columns=num_to_name)
        else:
            df = df.rename(columns=rename_dictionary)
    else:
        df = df.rename(columns=rename_dictionary)

    return df

def read_bim(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_bim")
    return df

def read_fam(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_fam")
    return df

def read_psam(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_psam")
    return df

def read_pvar(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_pvar")
    return df

def read_bgen_sample(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="bgen_sample")
    return df