#!/usr/bin/env python3
"""
Example script demonstrating how to use plot_arc to visualize BEDPE contact pairs.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from gwaslab.viz.viz_plot_arc import plot_arc

# Example 1: Plot arcs for a specific region
print("Example 1: Plotting arcs for chr10:127000000-128000000")
ax, df = plot_arc(
    'test_sim.bedpe.gz',
    region=('chr10', 127000000, 128000000),
    max_arcs=100,
    color="#FF0000",
    alpha=0.5,
    linewidth=1.0,
    arc_height=0.2,
    verbose=True
)
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('Track')
ax.set_title('BEDPE Contact Arcs - chr10:127-128Mb')
plt.savefig('example_arc_region.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved example_arc_region.png with {len(df)} arcs\n")

# Example 2: Plot all arcs (limited to 50 for visualization)
print("Example 2: Plotting first 50 arcs from entire file")
ax, df = plot_arc(
    'test_sim.bedpe.gz',
    max_arcs=50,
    color="#020080",
    alpha=0.3,
    linewidth=0.5,
    arc_height=0.1,
    verbose=True
)
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('Track')
ax.set_title('BEDPE Contact Arcs - All Pairs (sample)')
plt.savefig('example_arc_all.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved example_arc_all.png with {len(df)} arcs\n")

# Example 3: Color arcs by score (if score column exists)
print("Example 3: Plotting arcs colored by score")
ax, df = plot_arc(
    'test_sim.bedpe.gz',
    region=('chr10', 127000000, 128000000),
    max_arcs=100,
    score_col='score',
    color_by_score=True,
    cmap='plasma',
    alpha=0.6,
    linewidth=1.0,
    arc_height=0.2,
    verbose=True
)
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('Track')
ax.set_title('BEDPE Contact Arcs - Colored by Score')
plt.savefig('example_arc_score.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved example_arc_score.png with {len(df)} arcs\n")

print("All examples completed!")
