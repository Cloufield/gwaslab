"""
Materialize deferred AlphaGenome panels (ag_spec) via gwaslab-alphagenome.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_panel import Panel


def _locus_key(panel: Panel) -> Tuple[Any, ...]:
    vc = panel.get_kwarg("variant_context")
    vc_id = id(vc) if vc is not None else None
    region = panel.get_kwarg("region")
    build = panel.get_kwarg("build", "38")
    return (tuple(region) if region else None, str(build), vc_id)


def materialize_ag_panels(
    panels: List[Panel],
    log: Log = Log(),
    verbose: bool = True,
) -> List[Panel]:
    """
    Resolve ``ag_spec`` on panels by calling ``gwaslab_alphagenome.extract_batch``.

    Panels without ``ag_spec`` or with ``bundle`` already set are unchanged.
    """
    deferred: List[Tuple[int, Panel]] = []
    for i, p in enumerate(panels):
        if p.get_kwarg("ag_spec") is not None and p.get_kwarg("bundle") is None:
            deferred.append((i, p))

    if not deferred:
        return panels

    try:
        import gwaslab_alphagenome as glag
    except ImportError as e:
        raise ImportError(
            "Panels use ag_spec but gwaslab-alphagenome is not installed. "
            "Install with: pip install gwaslab-alphagenome"
        ) from e

    groups: Dict[Tuple[Any, ...], List[Tuple[int, Panel]]] = {}
    for idx, p in deferred:
        key = _locus_key(p)
        groups.setdefault(key, []).append((idx, p))

    out = list(panels)
    for _key, group in groups.items():
        ref_panel = group[0][1]
        region = ref_panel.get_kwarg("region")
        if region is None:
            raise ValueError("Deferred ag_spec panels require 'region' in panel kwargs")
        specs = [p.get_kwarg("ag_spec") for _, p in group]
        variant_context = ref_panel.get_kwarg("variant_context")
        build = ref_panel.get_kwarg("build", "38")
        log.write(
            f" -Materializing {len(specs)} AlphaGenome panel(s) via extract_batch...",
            verbose=verbose,
        )
        bundles = glag.extract_batch(
            region=region,
            specs=specs,
            variant_context=variant_context,
            build=build,
        )
        if len(bundles) != len(group):
            raise RuntimeError(
                f"extract_batch returned {len(bundles)} bundles for {len(group)} specs"
            )
        for (idx, panel), bundle in zip(group, bundles):
            panel.set_kwarg("bundle", bundle)
            panel.kwargs.pop("ag_spec", None)
            out[idx] = panel
    return out
