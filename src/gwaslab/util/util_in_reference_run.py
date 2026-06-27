"""
Reference-run planning: storage/CPU profiles, time estimates, progress tracking.

Baseline heuristics assume local SSD + baseline-tier 8-core desktop.
Calibration JSON under ~/.gwaslab/calibration/runs/, repo benchmarks/, or bundled
package data overrides heuristics when available.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

# --- Tunable SSD + baseline CPU constants (seconds) ---
_BCF_LOOKUP_SEC_PER_1K_TARGETS = 0.045
_BCF_LOOKUP_SEC_PER_GB_REF_UNINDEXED = 8.0
_BCF_LOOKUP_PER_CHR_BASE_SEC = 3.0
_BCF_LOOKUP_CONVERT_BCF_SEC_PER_GB = 45.0
_LOOKUP_ASSIGN_SEC_PER_1M_ROWS = 4.0
_TABIX_SEC_PER_VARIANT = 0.0025
_TSV_SCAN_SEC_PER_GB = 90.0
_FASTA_LOAD_SEC_PER_MB = 0.04
_FASTA_CHECK_SEC_PER_1K_VARIANTS = 0.002
_CALIBRATION_SCALE_MAX = 50.0
_TARGET_LIST_OPS = frozenset({"bcf_lookup_extract", "tabix_vcf_query"})

_IO_WEIGHTS: Dict[str, Tuple[float, float]] = {
    "bcf_lookup_extract": (0.6, 0.4),
    "lookup_assign": (0.25, 0.75),
    "tabix_vcf_query": (0.5, 0.5),
    "tsv_chunked_scan": (0.85, 0.15),
    "fasta_load_check": (0.55, 0.45),
    "harmonize_reference": (0.65, 0.35),
}

_STORAGE_IO_MULT: Dict[str, Tuple[float, float]] = {
    "ssd": (1.0, 1.0),
    "hdd": (2.5, 3.5),
    "network": (4.0, 8.0),
    "unknown": (1.5, 4.0),
}

_CPU_TIER_MULT: Dict[str, Tuple[float, float]] = {
    "fast": (0.7, 0.85),
    "baseline": (1.0, 1.0),
    "slow": (1.5, 2.5),
    "unknown": (1.0, 1.8),
}

_NETWORK_MOUNT_TYPES = frozenset({
    "nfs", "nfs4", "cifs", "fuse.sshfs", "fuse.glusterfs", "fuse.rclone",
})

_CALIBRATION_CACHE: Optional[List[dict]] = None
_COMPUTE_PROFILE_CACHE: Optional["ComputeProfile"] = None


class StorageProfile(str, Enum):
    SSD = "ssd"
    HDD = "hdd"
    NETWORK = "network"
    UNKNOWN = "unknown"


class CpuTier(str, Enum):
    SLOW = "slow"
    BASELINE = "baseline"
    FAST = "fast"
    UNKNOWN = "unknown"


@dataclass
class ComputeProfile:
    logical_cores: int
    physical_cores: Optional[int] = None
    threads_per_core: Optional[int] = None
    cpu_model: Optional[str] = None
    cpu_mhz: Optional[float] = None
    tier: CpuTier = CpuTier.BASELINE
    tier_override: Optional[str] = None
    ram_total_gb: Optional[float] = None
    ram_available_gb: Optional[float] = None

    def effective_workers(self, threads: int, task_count: int) -> int:
        t = max(1, int(threads))
        cores = max(1, self.logical_cores)
        tasks = max(1, int(task_count))
        return max(1, min(t, cores, tasks))

    def tier_label(self) -> str:
        if self.tier_override:
            return f"{self.tier.value} (user override)"
        return f"{self.tier.value} (detected)"


@dataclass
class ReferenceFileInfo:
    path: str
    kind: str = "unknown"
    size_bytes: int = 0
    indexed: bool = False
    record_count: Optional[int] = None
    storage_profile: StorageProfile = StorageProfile.UNKNOWN
    storage_override: Optional[str] = None
    device: Optional[str] = None
    rotational: Optional[int] = None
    mount_fstype: Optional[str] = None

    def profile_label(self) -> str:
        if self.storage_override:
            return f"{self.storage_profile.value} (user override)"
        return f"{self.storage_profile.value} (detected)"


@dataclass
class RunStep:
    name: str
    detail: str = ""
    estimate_seconds_low: float = 0.0
    estimate_seconds_high: float = 0.0
    io_weight: float = 0.5
    cpu_weight: float = 0.5
    calibrated: bool = False


@dataclass
class ReferenceRunPlan:
    operation: str
    steps: List[RunStep] = field(default_factory=list)
    total_estimate_low: float = 0.0
    total_estimate_high: float = 0.0
    calibrated_from_n: int = 0
    compute_profile: Optional[ComputeProfile] = None
    reference_files: List[ReferenceFileInfo] = field(default_factory=list)

    def recompute_totals(self) -> None:
        self.total_estimate_low = sum(s.estimate_seconds_low for s in self.steps)
        self.total_estimate_high = sum(s.estimate_seconds_high for s in self.steps)


def format_bytes(n: int) -> str:
    if n < 1024:
        return f"{n} B"
    if n < 1024 ** 2:
        return f"{n / 1024:.1f} KB"
    if n < 1024 ** 3:
        return f"{n / 1024 ** 2:.1f} MB"
    return f"{n / 1024 ** 3:.2f} GB"


def format_duration(seconds: float) -> str:
    if seconds < 0 or not (seconds < 1e12):
        return "?"
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        m = int(seconds // 60)
        s = int(seconds % 60)
        return f"{m}m {s}s" if s else f"{m}m"
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    return f"{h}h {m}m" if m else f"{h}h"


def _resolve_storage_profile_override() -> Optional[str]:
    v = os.environ.get("GWASLAB_STORAGE_PROFILE", "").strip().lower()
    if v in ("ssd", "hdd", "network"):
        return v
    return None


def _resolve_cpu_tier_override() -> Optional[str]:
    v = os.environ.get("GWASLAB_CPU_TIER", "").strip().lower()
    if v in ("slow", "baseline", "fast"):
        return v
    return None


def _parse_cpu_model_tier(model: str, logical_cores: int) -> CpuTier:
    m = model.lower()
    slow_hints = ("celeron", "atom", "pentium", "t2.small", "t2.micro", "core i3-2", "core i3-3")
    fast_hints = ("xeon gold", "xeon platinum", "epyc", "ryzen 9", "m2 pro", "m2 max", "m3 pro", "m3 max")
    if any(h in m for h in slow_hints) or logical_cores <= 4:
        if any(h in m for h in slow_hints) or logical_cores <= 2:
            return CpuTier.SLOW
    if logical_cores >= 16 or any(h in m for h in fast_hints):
        return CpuTier.FAST
    return CpuTier.BASELINE


def detect_compute_profile(cpu_tier: Optional[str] = None) -> ComputeProfile:
    global _COMPUTE_PROFILE_CACHE
    override = cpu_tier or _resolve_cpu_tier_override()
    if _COMPUTE_PROFILE_CACHE is not None and override is None:
        return _COMPUTE_PROFILE_CACHE

    logical = os.cpu_count() or 1
    try:
        logical = len(os.sched_getaffinity(0))
    except (AttributeError, NotImplementedError, OSError):
        pass

    physical = None
    threads_per_core = None
    cpu_model = None
    cpu_mhz = None

    if os.path.exists("/proc/cpuinfo"):
        try:
            with open("/proc/cpuinfo", "r") as f:
                text = f.read()
            for line in text.splitlines():
                if line.startswith("model name") and cpu_model is None:
                    cpu_model = line.split(":", 1)[1].strip()
                if line.startswith("cpu MHz") and cpu_mhz is None:
                    try:
                        cpu_mhz = float(line.split(":", 1)[1].strip())
                    except ValueError:
                        pass
                if line.startswith("cpu cores") and physical is None:
                    try:
                        physical = int(line.split(":", 1)[1].strip())
                    except ValueError:
                        pass
            if physical and logical:
                threads_per_core = max(1, logical // max(1, physical))
        except OSError:
            pass

    if override:
        tier = CpuTier(override)
        tier_override = override
    else:
        tier = _parse_cpu_model_tier(cpu_model or "", logical)
        tier_override = None

    ram_total_gb = None
    ram_available_gb = None
    if os.path.exists("/proc/meminfo"):
        try:
            mem: Dict[str, int] = {}
            with open("/proc/meminfo", "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        mem[parts[0].rstrip(":")] = int(parts[1])
            if "MemTotal" in mem:
                ram_total_gb = mem["MemTotal"] / 1024 / 1024
            if "MemAvailable" in mem:
                ram_available_gb = mem["MemAvailable"] / 1024 / 1024
        except OSError:
            pass

    profile = ComputeProfile(
        logical_cores=logical,
        physical_cores=physical,
        threads_per_core=threads_per_core,
        cpu_model=cpu_model,
        cpu_mhz=cpu_mhz,
        tier=tier,
        tier_override=tier_override,
        ram_total_gb=ram_total_gb,
        ram_available_gb=ram_available_gb,
    )
    if override is None:
        _COMPUTE_PROFILE_CACHE = profile
    return profile


def _mount_info_for_path(path: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Return (device, fstype, mount_point) for path."""
    path = os.path.abspath(path)
    best_mp = ""
    best_dev = None
    best_fstype = None
    if not os.path.exists("/proc/mounts"):
        return None, None, None
    try:
        with open("/proc/mounts", "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                dev, mp, fstype = parts[0], parts[1], parts[2]
                if path.startswith(mp) and len(mp) >= len(best_mp):
                    best_mp = mp
                    best_dev = dev
                    best_fstype = fstype
    except OSError:
        pass
    return best_dev, best_fstype, best_mp or None


def _block_rotational(device: str) -> Optional[int]:
    if not device:
        return None
    dev = device
    if dev.startswith("/dev/"):
        dev = dev[5:]
    # Handle /dev/sda1 -> sda
    base = re.sub(r"p?\d+$", "", dev)
    if base.startswith("mapper/"):
        return None
    rot_path = f"/sys/block/{base}/queue/rotational"
    try:
        with open(rot_path, "r") as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None


def detect_storage_profile(path: str, storage_profile: Optional[str] = None) -> StorageProfile:
    override = storage_profile or _resolve_storage_profile_override()
    if override:
        return StorageProfile(override)

    if not path or not os.path.exists(path):
        return StorageProfile.UNKNOWN

    device, fstype, _ = _mount_info_for_path(path)
    if fstype and fstype.lower() in _NETWORK_MOUNT_TYPES:
        return StorageProfile.NETWORK
    if fstype and fstype.lower().startswith("fuse."):
        return StorageProfile.NETWORK

    rot = _block_rotational(device or "")
    if rot == 0:
        return StorageProfile.SSD
    if rot == 1:
        return StorageProfile.HDD
    return StorageProfile.UNKNOWN


def _ref_kind(path: str) -> str:
    low = path.lower()
    if low.endswith(".bcf") or low.endswith(".bcf.csi"):
        return "bcf"
    if ".vcf" in low:
        return "vcf"
    if low.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fa.bgz", ".fasta.bgz")):
        return "fasta"
    if low.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz", ".lookup.txt.gz")):
        return "tsv"
    return "unknown"


def _is_indexed(path: str) -> bool:
    return (
        os.path.exists(path + ".tbi")
        or os.path.exists(path + ".csi")
        or os.path.exists(path + ".bcf.csi")
    )


def _bcftools_record_count(path: str) -> Optional[int]:
    try:
        out = subprocess.check_output(
            ["bcftools", "index", "-n", path],
            stderr=subprocess.DEVNULL,
            timeout=30,
            text=True,
        )
        return int(out.strip())
    except (subprocess.SubprocessError, ValueError, FileNotFoundError, OSError):
        return None


def collect_reference_file_info(
    path: str,
    storage_profile: Optional[str] = None,
) -> ReferenceFileInfo:
    override = storage_profile or _resolve_storage_profile_override()
    size = os.path.getsize(path) if path and os.path.exists(path) else 0
    device, fstype, _ = _mount_info_for_path(path) if path else (None, None, None)
    rot = _block_rotational(device or "") if device else None
    profile = detect_storage_profile(path, storage_profile=override)
    kind = _ref_kind(path) if path else "unknown"
    indexed = _is_indexed(path) if path else False
    record_count = None
    if path and kind in ("vcf", "bcf") and indexed:
        record_count = _bcftools_record_count(path)

    return ReferenceFileInfo(
        path=path,
        kind=kind,
        size_bytes=size,
        indexed=indexed,
        record_count=record_count,
        storage_profile=profile,
        storage_override=override,
        device=device,
        rotational=rot,
        mount_fstype=fstype,
    )


def _hybrid_estimate(
    io_base: float,
    cpu_base: float,
    io_weight: float,
    cpu_weight: float,
    storage: StorageProfile,
    compute: ComputeProfile,
    effective_workers: int,
) -> Tuple[float, float]:
    io_lo, io_hi = _STORAGE_IO_MULT.get(storage.value, _STORAGE_IO_MULT["unknown"])
    cpu_lo, cpu_hi = _CPU_TIER_MULT.get(compute.tier.value, _CPU_TIER_MULT["unknown"])
    par = max(1, effective_workers)
    io_part_lo = io_base * io_weight * io_lo
    io_part_hi = io_base * io_weight * io_hi
    cpu_part_lo = cpu_base * cpu_weight * cpu_lo / par
    cpu_part_hi = cpu_base * cpu_weight * cpu_hi / par
    return io_part_lo + cpu_part_lo, io_part_hi + cpu_part_hi


def _user_calibration_runs() -> Path:
    from gwaslab.bd.bd_download import get_default_directory

    return Path(get_default_directory()) / "calibration" / "runs"


def _repo_calibration_runs() -> Optional[Path]:
    here = Path(__file__).resolve()
    for parent in [here.parents[3], here.parents[2], Path.cwd()]:
        candidate = parent / "benchmarks" / "reference_run_calibration" / "runs"
        if candidate.is_dir():
            return candidate
    return None


def _bundled_calibration_runs() -> Optional[Path]:
    bundled = Path(__file__).resolve().parents[1] / "data" / "calibration_runs"
    if (bundled / "ci").is_dir():
        return bundled
    return None


def _calibration_roots() -> List[Path]:
    env = os.environ.get("GWASLAB_CALIBRATION_ROOT")
    if env:
        return [Path(env)]
    roots: List[Path] = []
    user = _user_calibration_runs()
    if user.is_dir():
        roots.append(user)
    repo = _repo_calibration_runs()
    if repo is not None:
        roots.append(repo)
    bundled = _bundled_calibration_runs()
    if bundled is not None:
        roots.append(bundled)
    return roots


def _calibration_root() -> Path:
    """Default runs directory for writes (benchmark harness); may not exist yet."""
    env = os.environ.get("GWASLAB_CALIBRATION_ROOT")
    if env:
        return Path(env)
    return _user_calibration_runs()


def load_calibration_runs(force: bool = False) -> List[dict]:
    global _CALIBRATION_CACHE
    if _CALIBRATION_CACHE is not None and not force:
        return _CALIBRATION_CACHE
    records: List[dict] = []
    seen_run_ids: set = set()
    for root in _calibration_roots():
        if not root.is_dir():
            continue
        for pattern in ("**/*.json", "*.json"):
            for p in root.glob(pattern):
                try:
                    with open(p, "r") as f:
                        data = json.load(f)
                    if not isinstance(data, dict) or "operation" not in data:
                        continue
                    run_id = data.get("run_id")
                    if run_id:
                        if run_id in seen_run_ids:
                            continue
                        seen_run_ids.add(run_id)
                    records.append(data)
                except (json.JSONDecodeError, OSError):
                    continue
    _CALIBRATION_CACHE = records
    return records


def _parse_duration_token(token: str) -> Optional[float]:
    """Parse '9m 16s', '572.192s', or '45s' into seconds."""
    token = token.strip()
    if not token:
        return None
    total = 0.0
    for part in token.replace(",", "").split():
        if part.endswith("ms"):
            continue
        if part.endswith("h"):
            total += float(part[:-1]) * 3600
        elif part.endswith("m"):
            total += float(part[:-1]) * 60
        elif part.endswith("s"):
            total += float(part[:-1])
        else:
            try:
                total += float(part)
            except ValueError:
                return None
    return total if total > 0 else None


def _calibration_admits(
    run: dict,
    operation: str,
    target_variants_request: int,
) -> bool:
    wall = run.get("timing", {}).get("wall_seconds")
    if wall is None or wall < 5.0:
        return False
    wl = run.get("workload", {})
    tv = int(wl.get("target_variants") or wl.get("sumstats_rows") or 0)
    if tv <= 0:
        return False
    if operation in _TARGET_LIST_OPS:
        min_tv = min(10000, max(1, target_variants_request // 100))
        if tv < min_tv:
            return False
        lookup_rows = wl.get("lookup_rows")
        if lookup_rows is not None and lookup_rows == 0 and wall < 30.0:
            return False
    return True


def _calibrated_estimate(
    operation: str,
    ref_size: int,
    target_variants: int,
    storage: StorageProfile,
    compute: ComputeProfile,
    task_count: int = 1,
    threads: int = 1,
    ref_indexed: bool = True,
) -> Optional[Tuple[float, float, int]]:
    runs = load_calibration_runs()
    admitted: List[Tuple[dict, float, int, int]] = []
    for r in runs:
        if r.get("operation") != operation:
            continue
        st = r.get("storage", {})
        if st.get("storage_profile_detected") != storage.value:
            continue
        if not _calibration_admits(r, operation, target_variants):
            continue
        wl = r.get("workload", {})
        tv = int(wl.get("target_variants") or wl.get("sumstats_rows") or 0)
        rs = int(st.get("ref_size_bytes") or 0)
        wall = float(r.get("timing", {}).get("wall_seconds"))
        admitted.append((r, wall, rs, tv))

    if not admitted:
        return None

    def dist(entry: Tuple[dict, float, int, int]) -> float:
        _, _, rs, tv = entry
        d_tv = abs(tv - target_variants) / max(target_variants, tv, 1)
        if operation in _TARGET_LIST_OPS:
            return d_tv
        d_rs = abs(rs - ref_size) / max(ref_size, rs, 1)
        return d_rs + d_tv

    best_r, wall, rs, tv = min(admitted, key=dist)
    if tv > 0 and target_variants < tv / 5:
        return None
    wl = best_r.get("workload", {})
    threads_cal = int(
        wl.get("extract_threads_requested")
        or wl.get("threads_requested")
        or wl.get("threads_effective")
        or 1
    )
    task_count_cal = int(wl.get("task_count") or task_count or 1)
    eff_cal = compute.effective_workers(threads_cal, task_count_cal)
    eff_req = compute.effective_workers(threads, task_count)

    if operation in _TARGET_LIST_OPS:
        rate = wall / max(tv, 1)
        scaled = rate * target_variants * (eff_cal / max(eff_req, 1))
        if not ref_indexed and rs > 0 and ref_size > rs:
            scaled *= (ref_size / rs) ** 0.25
    else:
        scale = 1.0
        if rs > 0 and ref_size > 0:
            scale *= ref_size / rs
        if tv > 0 and target_variants > 0:
            scale *= target_variants / tv
        if rs > 0 and ref_size > 0 and tv > 0 and target_variants > 0:
            scale = ((ref_size / rs) * (target_variants / tv)) ** 0.5
        scaled = wall * scale

    if wall > 0 and scaled / wall > _CALIBRATION_SCALE_MAX:
        return None

    low = scaled * 0.75
    high = scaled * 1.5
    return low, high, len(admitted)


def parse_harmonize_log_steps(log_text: str) -> List[dict]:
    """Extract per-step wall times from a harmonize verbose log."""
    steps: List[dict] = []

    m_check = re.search(
        r"Start to check if NEA is aligned with reference sequence.*?Time taken: ([\d.]+)s",
        log_text,
        re.DOTALL,
    )
    if m_check:
        steps.append({"name": "check_ref", "wall_seconds": float(m_check.group(1))})

    extract_names = ("infer_extract", "rsid_extract")
    extract_idx = 0
    for m in re.finditer(
        r"Lookup extraction finished in ([^\(]+)\(([\d,]+) rows\)",
        log_text,
    ):
        dur = _parse_duration_token(m.group(1).strip()) or 0.0
        rows = int(m.group(2).replace(",", ""))
        name = extract_names[extract_idx] if extract_idx < len(extract_names) else f"extract_{extract_idx}"
        steps.append({"name": name, "wall_seconds": dur, "lookup_rows": rows})
        extract_idx += 1

    assign_names = ("infer_assign", "rsid_assign")
    assign_idx = 0
    for m in re.finditer(
        r"Finished assigning from lookup table\.\n.*?Time taken: ([\d.]+)s",
        log_text,
        re.DOTALL,
    ):
        name = assign_names[assign_idx] if assign_idx < len(assign_names) else f"assign_{assign_idx}"
        steps.append({"name": name, "wall_seconds": float(m.group(1))})
        assign_idx += 1

    return steps


def estimate_run_plan(
    operation: str,
    compute_profile: Optional[ComputeProfile] = None,
    storage_profile: Optional[StorageProfile] = None,
    ref_path: Optional[str] = None,
    ref_info: Optional[ReferenceFileInfo] = None,
    threads: int = 1,
    task_count: int = 1,
    target_variants: int = 0,
    sumstats_rows: int = 0,
    lookup_size_bytes: int = 0,
    chunksize: int = 5_000_000,
    convert_to_bcf: bool = False,
    chroms_loaded: int = 0,
    variants_to_check: int = 0,
    cpu_tier: Optional[str] = None,
    storage_profile_override: Optional[str] = None,
    log_subplans: bool = True,
) -> ReferenceRunPlan:
    compute = compute_profile or detect_compute_profile(cpu_tier=cpu_tier)
    info = ref_info
    if info is None and ref_path:
        info = collect_reference_file_info(ref_path, storage_profile=storage_profile_override)
    storage = storage_profile or (info.storage_profile if info else StorageProfile.UNKNOWN)
    ref_size = info.size_bytes if info else 0
    targets = target_variants or sumstats_rows
    eff = compute.effective_workers(threads, task_count)
    io_w, cpu_w = _IO_WEIGHTS.get(operation, (0.5, 0.5))

    plan = ReferenceRunPlan(
        operation=operation,
        compute_profile=compute,
        reference_files=[info] if info else [],
    )

    cal = _calibrated_estimate(
        operation,
        ref_size,
        targets,
        storage,
        compute,
        task_count=task_count,
        threads=threads,
        ref_indexed=getattr(info, "indexed", True) if info else True,
    )
    if cal:
        low, high, n = cal
        plan.steps.append(
            RunStep(
                name=operation,
                detail=f"calibrated from {n} local run(s)",
                estimate_seconds_low=low,
                estimate_seconds_high=high,
                io_weight=io_w,
                cpu_weight=cpu_w,
                calibrated=True,
            )
        )
        plan.calibrated_from_n = n
        plan.recompute_totals()
        return plan

    low, high = 0.0, 0.0
    detail_parts: List[str] = []

    if operation == "bcf_lookup_extract":
        io_base = targets / 1000 * _BCF_LOOKUP_SEC_PER_1K_TARGETS
        io_base += (task_count / max(eff, 1)) * _BCF_LOOKUP_PER_CHR_BASE_SEC
        if info is None or not getattr(info, "indexed", True):
            gb = ref_size / (1024 ** 3)
            io_base += gb * _BCF_LOOKUP_SEC_PER_GB_REF_UNINDEXED
        cpu_base = targets / 1000 * 0.02
        if convert_to_bcf:
            gb = ref_size / (1024 ** 3)
            io_base += gb * _BCF_LOOKUP_CONVERT_BCF_SEC_PER_GB
        detail_parts = [f"{task_count} chr tasks", f"{targets:,} targets", f"{eff}/{threads} workers"]
        low, high = _hybrid_estimate(io_base, cpu_base, io_w, cpu_w, storage, compute, eff)

    elif operation == "lookup_assign":
        lookup_rows_est = lookup_size_bytes // 80 if lookup_size_bytes else targets
        io_base = lookup_size_bytes / (1024 ** 3) * 30.0
        cpu_base = lookup_rows_est / 1e6 * _LOOKUP_ASSIGN_SEC_PER_1M_ROWS
        detail_parts = [f"~{lookup_rows_est:,} lookup rows est.", f"{sumstats_rows:,} sumstats rows"]
        low, high = _hybrid_estimate(io_base, cpu_base, io_w, cpu_w, storage, compute, 1)

    elif operation == "tabix_vcf_query":
        io_base = targets * _TABIX_SEC_PER_VARIANT * 0.6 + ref_size / (1024 ** 3) * 5.0
        cpu_base = targets * _TABIX_SEC_PER_VARIANT * 0.4
        detail_parts = [f"{targets:,} variants", f"{eff}/{threads} workers"]
        low, high = _hybrid_estimate(io_base, cpu_base, io_w, cpu_w, storage, compute, eff)

    elif operation == "tsv_chunked_scan":
        gb = ref_size / (1024 ** 3)
        blocks = max(1, int(gb * 1e9 / max(chunksize, 1)) + 1) if ref_size else 1
        io_base = gb * _TSV_SCAN_SEC_PER_GB
        cpu_base = blocks * 2.0
        detail_parts = [f"~{blocks} blocks", format_bytes(ref_size)]
        low, high = _hybrid_estimate(io_base, cpu_base, io_w, cpu_w, storage, compute, 1)

    elif operation == "fasta_load_check":
        mb = ref_size / (1024 ** 2)
        io_base = mb * _FASTA_LOAD_SEC_PER_MB * max(1, chroms_loaded) / max(1, chroms_loaded)
        cpu_base = variants_to_check / 1000 * _FASTA_CHECK_SEC_PER_1K_VARIANTS
        detail_parts = [f"{chroms_loaded} chroms", f"{variants_to_check:,} variants to check"]
        low, high = _hybrid_estimate(io_base, cpu_base, io_w, cpu_w, storage, compute, 1)

    elif operation == "harmonize_reference":
        # Composite — caller should pass pre-built steps or we sum sub-ops
        detail_parts = ["combined reference steps"]
        low, high = 60.0, 600.0

    else:
        detail_parts = ["unknown operation"]
        low, high = 30.0, 300.0

    plan.steps.append(
        RunStep(
            name=operation,
            detail=", ".join(detail_parts),
            estimate_seconds_low=max(0.5, low),
            estimate_seconds_high=max(low * 1.2, high),
            io_weight=io_w,
            cpu_weight=cpu_w,
        )
    )
    plan.recompute_totals()
    return plan


def plan_harmonize_reference_steps(
    ref_seq: Optional[str] = None,
    ref_infer: Optional[str] = None,
    ref_rsid_vcf: Optional[str] = None,
    ref_rsid_tsv: Optional[str] = None,
    sweep_mode: bool = False,
    threads: int = 1,
    extract_threads: Optional[int] = None,
    sumstats_rows: int = 0,
    target_infer: int = 0,
    target_rsid: int = 0,
    variants_check_ref: int = 0,
    chroms_in_sumstats: int = 0,
    cpu_tier: Optional[str] = None,
    storage_profile: Optional[str] = None,
) -> ReferenceRunPlan:
    compute = detect_compute_profile(cpu_tier=cpu_tier)
    combined = ReferenceRunPlan(operation="harmonize_reference", compute_profile=compute)
    extract_thr = extract_threads if extract_threads is not None else 1
    chr_tasks = max(1, chroms_in_sumstats)

    if ref_seq and variants_check_ref > 0:
        info = collect_reference_file_info(ref_seq, storage_profile=storage_profile)
        combined.reference_files.append(info)
        sub = estimate_run_plan(
            "fasta_load_check",
            compute_profile=compute,
            ref_info=info,
            chroms_loaded=max(1, chroms_in_sumstats),
            variants_to_check=variants_check_ref,
            storage_profile_override=storage_profile,
        )
        for s in sub.steps:
            s.name = "check_ref"
            combined.steps.append(s)

    if ref_infer and target_infer > 0:
        info = collect_reference_file_info(ref_infer, storage_profile=storage_profile)
        combined.reference_files.append(info)
        op = "bcf_lookup_extract" if sweep_mode else "tabix_vcf_query"
        sub = estimate_run_plan(
            op,
            compute_profile=compute,
            ref_info=info,
            threads=extract_thr if sweep_mode else threads,
            task_count=chr_tasks if sweep_mode else max(1, threads),
            target_variants=target_infer if sweep_mode else target_infer,
            storage_profile_override=storage_profile,
        )
        for s in sub.steps:
            s.name = "infer_strand" + (" (sweep)" if sweep_mode else " (tabix)")
            combined.steps.append(s)
        if sweep_mode:
            sub2 = estimate_run_plan(
                "lookup_assign",
                compute_profile=compute,
                target_variants=target_infer,
                sumstats_rows=sumstats_rows,
            )
            for s in sub2.steps:
                s.name = "infer_strand assign"
                combined.steps.append(s)

    if ref_rsid_tsv and target_rsid > 0:
        info = collect_reference_file_info(ref_rsid_tsv, storage_profile=storage_profile)
        combined.reference_files.append(info)
        sub = estimate_run_plan(
            "tsv_chunked_scan",
            compute_profile=compute,
            ref_info=info,
            target_variants=target_rsid,
            storage_profile_override=storage_profile,
        )
        for s in sub.steps:
            s.name = "assign_rsid (TSV)"
            combined.steps.append(s)

    if ref_rsid_vcf and target_rsid > 0:
        info = collect_reference_file_info(ref_rsid_vcf, storage_profile=storage_profile)
        combined.reference_files.append(info)
        if sweep_mode:
            sub = estimate_run_plan(
                "bcf_lookup_extract",
                compute_profile=compute,
                ref_info=info,
                threads=extract_thr,
                task_count=chr_tasks,
                target_variants=target_rsid,
                storage_profile_override=storage_profile,
            )
            for s in sub.steps:
                s.name = "assign_rsid (sweep extract)"
                combined.steps.append(s)
            sub2 = estimate_run_plan(
                "lookup_assign",
                compute_profile=compute,
                target_variants=target_rsid,
                sumstats_rows=sumstats_rows,
            )
            for s in sub2.steps:
                s.name = "assign_rsid (sweep assign)"
                combined.steps.append(s)
        else:
            sub = estimate_run_plan(
                "tabix_vcf_query",
                compute_profile=compute,
                ref_info=info,
                threads=threads,
                task_count=max(1, threads),
                target_variants=target_rsid,
                storage_profile_override=storage_profile,
            )
            for s in sub.steps:
                s.name = "assign_rsid (tabix)"
                combined.steps.append(s)

    combined.calibrated_from_n = sum(1 for s in combined.steps if s.calibrated)
    combined.recompute_totals()
    return combined


class RunProgressTracker:
    def __init__(self, total: int, label: str = "items") -> None:
        self.total = max(1, total)
        self.label = label
        self.completed = 0
        self._start = time.time()
        self._last_label: Optional[str] = None

    @property
    def elapsed(self) -> float:
        return time.time() - self._start

    def update(self, label: Optional[str] = None) -> None:
        self.completed += 1
        if label is not None:
            self._last_label = label

    def eta_seconds(self) -> Optional[float]:
        if self.completed <= 0:
            return None
        rate = self.completed / self.elapsed
        remaining = self.total - self.completed
        if rate <= 0:
            return None
        return remaining / rate

    def log_progress(self, log: Any, verbose: bool = True) -> None:
        if hasattr(log, "log_progress"):
            log.log_progress(
                done=self.completed,
                total=self.total,
                label=self._last_label or self.label,
                elapsed_s=self.elapsed,
                eta_s=self.eta_seconds(),
                verbose=verbose,
            )


class MemoryTracker:
    """Track peak RSS during a block."""

    def __init__(self) -> None:
        self._start_rss = self._rss_mb()
        self.peak_rss_mb = self._start_rss

    @staticmethod
    def _rss_mb() -> float:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        except Exception:
            return 0.0

    def sample(self) -> None:
        rss = self._rss_mb()
        if rss > self.peak_rss_mb:
            self.peak_rss_mb = rss

    def delta_mb(self) -> float:
        return max(0.0, self.peak_rss_mb - self._start_rss)


def snapshot_machine_for_json() -> dict:
    compute = detect_compute_profile()
    return {
        "hostname": os.environ.get("HOSTNAME", ""),
        "cpu_model": compute.cpu_model,
        "logical_cores": compute.logical_cores,
        "physical_cores": compute.physical_cores,
        "threads_per_core": compute.threads_per_core,
        "cpu_mhz_current": compute.cpu_mhz,
        "cpu_tier_detected": compute.tier.value,
        "cpu_tier_override": compute.tier_override,
        "ram_total_gb": compute.ram_total_gb,
        "ram_available_gb_at_start": compute.ram_available_gb,
    }


def snapshot_storage_for_json(path: str, storage_override: Optional[str] = None) -> dict:
    info = collect_reference_file_info(path, storage_profile=storage_override)
    return {
        "ref_path": path,
        "ref_size_bytes": info.size_bytes,
        "ref_kind": info.kind,
        "indexed": info.indexed,
        "record_count": info.record_count,
        "workspace_device": info.device,
        "rotational": info.rotational,
        "storage_profile_detected": info.storage_profile.value,
        "storage_profile_override": info.storage_override,
        "mount_fstype": info.mount_fstype,
    }
