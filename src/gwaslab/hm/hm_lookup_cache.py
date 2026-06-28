"""Persistent Parquet lookup cache for VCF/BCF extract (memory-safe, per-chr partitions).
"""
from __future__ import annotations

import json
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Set, Tuple

import pandas as pd

LOOKUP_CACHE_FORMAT = "parquet_chr_partition"
SCHEMA_VERSION = 1


def is_lookup_bundle(path: str | Path) -> bool:
    p = Path(path)
    return p.is_dir() and (p / "meta.json").is_file()


def _ref_stem(vcf_path: str) -> str:
    name = Path(vcf_path).name
    for suffix in (".vcf.gz", ".bcf", ".vcf", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def lookup_cache_bundle_path(vcf_path: str, assign_cols: Sequence[str]) -> Path:
    """Cache root: ~/.gwaslab/lookup/{ref_stem}.lookup.{cols_tag}/
"""
    from gwaslab.bd.bd_download import get_default_directory

    cols_tag = "-".join(sorted(assign_cols))
    base = Path(get_default_directory()) / "lookup"
    base.mkdir(parents=True, exist_ok=True)
    return base / f"{_ref_stem(vcf_path)}.lookup.{cols_tag}"


def _chr_tag(chr_val) -> str:
    return f"chr{int(chr_val)}"


def _targets_path(bundle: Path, chr_val) -> Path:
    return bundle / "targets" / f"{_chr_tag(chr_val)}.parquet"


def _data_base_path(bundle: Path, chr_val) -> Path:
    return bundle / "data" / f"{_chr_tag(chr_val)}.parquet"


def _data_delta_paths(bundle: Path, chr_val) -> list[Path]:
    data_dir = bundle / "data"
    prefix = f"{_chr_tag(chr_val)}.delta."
    return sorted(data_dir.glob(f"{prefix}*.parquet")) if data_dir.is_dir() else []


def _read_cached_pos_set(bundle: Path, chr_val) -> Set[int]:
    path = _targets_path(bundle, chr_val)
    if not path.is_file():
        return set()
    pos = pd.read_parquet(path, columns=["POS"])["POS"]
    return set(int(x) for x in pos.dropna().unique())


def _write_pos_parquet(path: Path, positions: Set[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not positions:
        if path.is_file():
            path.unlink()
        return
    df = pd.DataFrame({"POS": sorted(positions)})
    df.to_parquet(path, compression="zstd", index=False)


def lookup_cache_meta_valid(
    bundle: Path,
    vcf_path: str,
    assign_cols: Sequence[str],
) -> bool:
    meta_path = bundle / "meta.json"
    if not meta_path.is_file():
        return False
    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except (json.JSONDecodeError, OSError):
        return False
    ref = Path(vcf_path).resolve()
    if not ref.is_file():
        return False
    try:
        st = ref.stat()
    except OSError:
        return False
    if meta.get("schema_version") != SCHEMA_VERSION:
        return False
    if meta.get("lookup_cache_format") != LOOKUP_CACHE_FORMAT:
        return False
    if list(meta.get("assign_cols", [])) != list(assign_cols):
        return False
    if meta.get("ref_path") != str(ref):
        return False
    if meta.get("ref_size") is not None and meta.get("ref_size") != st.st_size:
        return False
    if meta.get("ref_mtime") is not None and meta.get("ref_mtime") != int(st.st_mtime):
        return False
    return True


def _write_meta(bundle: Path, vcf_path: str, assign_cols: Sequence[str], lookup_rows: int) -> None:
    ref = Path(vcf_path).resolve()
    ref_size: Optional[int] = None
    ref_mtime: Optional[int] = None
    if ref.is_file():
        st = ref.stat()
        ref_size = st.st_size
        ref_mtime = int(st.st_mtime)
    meta = {
        "schema_version": SCHEMA_VERSION,
        "lookup_cache_format": LOOKUP_CACHE_FORMAT,
        "ref_path": str(ref),
        "ref_size": ref_size,
        "ref_mtime": ref_mtime,
        "assign_cols": list(assign_cols),
        "lookup_rows": lookup_rows,
        "updated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    bundle.mkdir(parents=True, exist_ok=True)
    with open(bundle / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)


def _invalidate_bundle(bundle: Path) -> None:
    if bundle.is_dir():
        shutil.rmtree(bundle)


def compute_target_delta_by_chr(
    targets: pd.DataFrame,
    bundle: Path,
) -> pd.DataFrame:
    """Return CHR,POS rows not yet in bundle targets/ (per-chr POS sets).
"""
    if targets.empty:
        return targets.copy()
    work = targets[["CHR", "POS"]].drop_duplicates()
    work["CHR"] = pd.to_numeric(work["CHR"], errors="coerce")
    work["POS"] = pd.to_numeric(work["POS"], errors="coerce")
    work = work.dropna(subset=["CHR", "POS"])
    delta_parts: list[pd.DataFrame] = []
    for chr_val, grp in work.groupby("CHR", sort=True):
        c = int(chr_val)
        current = {int(x) for x in grp["POS"].unique()}
        cached = _read_cached_pos_set(bundle, c)
        missing = current - cached
        if missing:
            delta_parts.append(
                pd.DataFrame({"CHR": c, "POS": sorted(missing)})
            )
    if not delta_parts:
        return pd.DataFrame(columns=["CHR", "POS"])
    return pd.concat(delta_parts, ignore_index=True)


def init_bundle_from_extract_build(
    build_path: str,
    bundle: Path,
    vcf_path: str,
    assign_cols: Sequence[str],
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    _invalidate_bundle(bundle)
    (bundle / "data").mkdir(parents=True, exist_ok=True)
    (bundle / "targets").mkdir(parents=True, exist_ok=True)

    writers: dict[int, pq.ParquetWriter] = {}
    pos_sets: dict[int, Set[int]] = {}
    total_rows = 0
    chunksize = 1_000_000

    for chunk in pd.read_csv(build_path, sep="\t", chunksize=chunksize):
        if chunk.empty:
            continue
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["POS"] = pd.to_numeric(chunk["POS"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "POS"])
        for chr_val, grp in chunk.groupby("CHR", sort=True):
            c = int(chr_val)
            path = _data_base_path(bundle, c)
            table = pa.Table.from_pandas(grp, preserve_index=False)
            if c not in writers:
                writers[c] = pq.ParquetWriter(
                    str(path), table.schema, compression="zstd"
                )
            writers[c].write_table(table)
            total_rows += len(grp)
            s = pos_sets.setdefault(c, set())
            s.update(int(x) for x in grp["POS"].unique())

    for w in writers.values():
        w.close()

    for c, positions in pos_sets.items():
        _write_pos_parquet(_targets_path(bundle, c), positions)

    _write_meta(bundle, vcf_path, assign_cols, total_rows)


def append_delta_to_bundle(
    bundle: Path,
    delta_parquet: str,
    vcf_path: str,
    assign_cols: Sequence[str],
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / "data").mkdir(parents=True, exist_ok=True)
    (bundle / "targets").mkdir(parents=True, exist_ok=True)

    delta_id = datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    writers: dict[int, pq.ParquetWriter] = {}
    new_pos: dict[int, Set[int]] = {}
    added_rows = 0
    chunksize = 500_000

    def _ingest_chunk(chunk: pd.DataFrame) -> None:
        nonlocal added_rows
        if chunk.empty:
            return
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["POS"] = pd.to_numeric(chunk["POS"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "POS"])
        for chr_val, grp in chunk.groupby("CHR", sort=True):
            c = int(chr_val)
            delta_path = bundle / "data" / f"{_chr_tag(c)}.delta.{delta_id}.parquet"
            table = pa.Table.from_pandas(grp, preserve_index=False)
            if c not in writers:
                writers[c] = pq.ParquetWriter(
                    str(delta_path), table.schema, compression="zstd"
                )
            writers[c].write_table(table)
            added_rows += len(grp)
            s = new_pos.setdefault(c, set())
            s.update(int(x) for x in grp["POS"].unique())

    is_parquet = str(delta_parquet).lower().endswith(".parquet")
    if is_parquet:
        try:
            pf = pq.ParquetFile(delta_parquet)
            for batch in pf.iter_batches(batch_size=chunksize):
                _ingest_chunk(batch.to_pandas())
        except Exception:
            is_parquet = False
    if not is_parquet:
        for chunk in pd.read_csv(delta_parquet, sep="\t", chunksize=chunksize):
            _ingest_chunk(chunk)

    for w in writers.values():
        w.close()

    for c, positions in new_pos.items():
        cached = _read_cached_pos_set(bundle, c)
        cached.update(positions)
        _write_pos_parquet(_targets_path(bundle, c), cached)

    meta_path = bundle / "meta.json"
    prev_rows = 0
    if meta_path.is_file():
        try:
            with open(meta_path) as f:
                prev_rows = int(json.load(f).get("lookup_rows", 0))
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            prev_rows = 0
    _write_meta(bundle, vcf_path, assign_cols, prev_rows + added_rows)


def lookup_bundle_size_bytes(bundle: Path) -> int:
    total = 0
    if not bundle.is_dir():
        return 0
    for root, _, files in os.walk(bundle):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return total


def iter_lookup_bundle_chunks(
    bundle: Path,
    usecols: Sequence[str],
    chrs: Set[int],
    pos_filter_by_chr: Optional[Dict[int, Set[int]]],
    dtype: dict,
    chunksize: int,
) -> Iterable[pd.DataFrame]:
    import pyarrow.parquet as pq

    usecols_list = list(usecols)
    for c in sorted(chrs):
        paths: list[Path] = []
        base = _data_base_path(bundle, c)
        if base.is_file():
            paths.append(base)
        paths.extend(_data_delta_paths(bundle, c))
        pos_allow = pos_filter_by_chr.get(c) if pos_filter_by_chr else None
        for path in paths:
            pf = pq.ParquetFile(str(path))
            cols = [x for x in usecols_list if x in pf.schema.names]
            if not cols:
                continue
            for batch in pf.iter_batches(batch_size=chunksize, columns=cols):
                chunk = batch.to_pandas()
                if pos_allow is not None:
                    chunk["POS"] = pd.to_numeric(chunk["POS"], errors="coerce")
                    chunk = chunk[chunk["POS"].isin(pos_allow)]
                    if chunk.empty:
                        continue
                for col, col_dtype in dtype.items():
                    if col not in chunk.columns:
                        continue
                    if col_dtype == "category":
                        chunk[col] = chunk[col].astype("category")
                    else:
                        chunk[col] = chunk[col].astype(col_dtype)
                yield chunk


def lookup_bundle_columns(bundle: Path) -> list[str]:
    import pyarrow.parquet as pq

    data_dir = bundle / "data"
    if not data_dir.is_dir():
        return []
    for path in sorted(data_dir.glob("*.parquet")):
        if ".delta." in path.name:
            continue
        return list(pq.ParquetFile(str(path)).schema.names)
    for path in sorted(data_dir.glob("*.parquet")):
        return list(pq.ParquetFile(str(path)).schema.names)
    return []


def lookup_bundle_has_rows(bundle: Path) -> bool:
    data_dir = bundle / "data"
    if not data_dir.is_dir():
        return False
    for path in data_dir.glob("*.parquet"):
        try:
            import pyarrow.parquet as pq

            if pq.ParquetFile(str(path)).metadata.num_rows > 0:
                return True
        except Exception:
            continue
    return False


def resolve_lookup_for_vcf(
    vcf_path: str,
    targets: pd.DataFrame,
    assign_cols: Sequence[str],
    mapper,
    lookup_path: Optional[str] = None,
    tsv_path: Optional[str] = None,
    reuse_lookup: bool = True,
    threads: int = 6,
    extract_threads: Optional[int] = None,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    verbose: bool = True,
    log_run_plan: bool = True,
    log=None,
) -> Tuple[str, bool]:
    """Resolve lookup table path; populate cache under ~/.gwaslab/lookup when needed.
"""
    import os
    import tempfile

    from gwaslab.hm.hm_assign_rsid import (
        _convert_vcf_to_bcf,
        _extract_lookup_table_from_vcf_bcf,
        is_vcf_file,
    )
    from gwaslab.info.g_Log import Log

    if log is None:
        log = Log()

    assign_cols_list = list(assign_cols)

    if lookup_path and os.path.exists(lookup_path):
        log.write(f" -Using lookup_path: {lookup_path}", verbose=verbose)
        return lookup_path, False

    if tsv_path and reuse_lookup and os.path.exists(tsv_path):
        log.write(f" -Reusing lookup TSV: {tsv_path}", verbose=verbose)
        return tsv_path, False

    path_to_use = vcf_path
    if convert_to_bcf and is_vcf_file(path_to_use):
        if len(path_to_use) < 4 or path_to_use[-4:] != ".bcf":
            log.write(
                f" -Converting VCF to BCF (strip_info={strip_info})...",
                verbose=verbose,
            )
            path_to_use = _convert_vcf_to_bcf(
                path_to_use, threads=threads, strip=strip_info, log=log, verbose=verbose
            )

    mapper.detect_reference_format(path_to_use)
    # Keep sumstats CHR notation; _extract_lookup_table_from_vcf_bcf converts to reference.
    work_targets = targets[["CHR", "POS"]].copy()
    work_targets = work_targets.dropna(subset=["CHR", "POS"]).drop_duplicates()

    bundle = lookup_cache_bundle_path(path_to_use, assign_cols_list)

    if reuse_lookup and is_lookup_bundle(bundle) and lookup_cache_meta_valid(
        bundle, path_to_use, assign_cols_list
    ):
        delta = compute_target_delta_by_chr(work_targets, bundle)
        if delta.empty:
            log.write(
                f" -Lookup cache hit ({bundle}), skipping bcftools extract",
                verbose=verbose,
            )
            return str(bundle), False
        log.write(
            f" -Lookup cache partial hit: {len(delta):,} new CHR:POS to extract",
            verbose=verbose,
        )
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".lookup.parquet")
        tmp_path = tmp.name
        tmp.close()
        try:
            _extract_lookup_table_from_vcf_bcf(
                vcf_path=path_to_use,
                sumstats=delta,
                mapper=mapper,
                assign_cols=assign_cols_list,
                out_lookup=tmp_path,
                threads=threads,
                extract_threads=extract_threads,
                verbose=verbose,
                log=log,
                log_run_plan=log_run_plan,
            )
            append_delta_to_bundle(bundle, tmp_path, path_to_use, assign_cols_list)
        finally:
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        return str(bundle), False

    if reuse_lookup and is_lookup_bundle(bundle):
        log.write(
            " -Lookup cache meta invalid (ref changed), rebuilding bundle...",
            verbose=verbose,
        )
        _invalidate_bundle(bundle)

    log.write(
        f" -Building lookup cache at {bundle} (full extract)...",
        verbose=verbose,
    )
    build_tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".lookup.build.tsv")
    build_path = build_tmp.name
    build_tmp.close()
    try:
        _extract_lookup_table_from_vcf_bcf(
            vcf_path=path_to_use,
            sumstats=work_targets,
            mapper=mapper,
            assign_cols=assign_cols_list,
            out_lookup=build_path,
            threads=threads,
            extract_threads=extract_threads,
            verbose=verbose,
            log=log,
            log_run_plan=log_run_plan,
        )
        init_bundle_from_extract_build(
            build_path, bundle, path_to_use, assign_cols_list
        )
    finally:
        try:
            os.remove(build_path)
        except OSError:
            pass
    return str(bundle), False
