"""
GWASLab-native track bundle dataclasses for AlphaGenome (and other) predicted tracks.

The gwaslab-alphagenome wrapper converts API output into these types; GWASLab plots them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

RegionTuple = Tuple[int, int, int]


@dataclass
class TrackBundle:
    """1D predicted signal track(s)."""

    values: np.ndarray
    metadata: pd.DataFrame
    resolution: int
    region: RegionTuple
    source: str = "alphagenome"

    def __post_init__(self) -> None:
        if self.values.ndim == 1:
            self.values = self.values.reshape(-1, 1)
        if "name" not in self.metadata.columns or "strand" not in self.metadata.columns:
            raise ValueError("TrackBundle metadata must include 'name' and 'strand' columns")

    @property
    def num_tracks(self) -> int:
        return int(self.values.shape[1])

    @property
    def num_axes(self) -> int:
        return self.num_tracks


@dataclass
class OverlayBundle:
    """REF/ALT paired tracks from variant effect prediction."""

    tracks: Dict[str, TrackBundle]
    region: RegionTuple
    variant_pos: Optional[int] = None
    variant_label: Optional[str] = None

    @property
    def num_axes(self) -> int:
        if not self.tracks:
            return 1
        first = next(iter(self.tracks.values()))
        return first.num_axes


@dataclass
class ContactBundle:
    """2D contact map track(s)."""

    values: np.ndarray
    metadata: pd.DataFrame
    resolution: int
    region: RegionTuple
    source: str = "alphagenome"

    def __post_init__(self) -> None:
        if self.values.ndim == 2:
            self.values = self.values[:, :, np.newaxis]
        if "name" not in self.metadata.columns:
            raise ValueError("ContactBundle metadata must include a 'name' column")

    @property
    def num_tracks(self) -> int:
        return int(self.values.shape[2])

    @property
    def num_axes(self) -> int:
        return self.num_tracks


@dataclass
class JunctionBundle:
    """Splice junction sashimi data."""

    junctions: pd.DataFrame
    values: np.ndarray
    metadata: pd.DataFrame
    region: RegionTuple
    strands: Tuple[str, ...] = ("+", "-")
    source: str = "alphagenome"

    @property
    def num_tracks(self) -> int:
        return len(self.metadata)

    @property
    def num_axes(self) -> int:
        return self.num_tracks * len(self.strands)


Bundle = Union[TrackBundle, OverlayBundle, ContactBundle, JunctionBundle]


def bundle_num_axes(bundle: Bundle) -> int:
    return bundle.num_axes


def x_positions(bundle: TrackBundle) -> np.ndarray:
    """Genomic x coordinates (bp) for each bin center."""
    start = bundle.region[1]
    n = bundle.values.shape[0]
    return np.arange(n) * bundle.resolution + start + bundle.resolution / 2.0
