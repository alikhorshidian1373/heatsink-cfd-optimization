"""
Shared plotting style for the heat-sink CFD study.

One place to change fonts, colours and figure geometry so every figure in the
report reads as part of the same system.

Palette slots are a validated categorical set (colour-blind separation checked):
    slot 1  blue    #2a78d6
    slot 2  orange  #eb6834
    slot 3  aqua    #1baf7a
Assign them in fixed order; never cycle to a generated hue.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- palette -----------------------------------------------------------------
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
SERIES = (BLUE, ORANGE, AQUA)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8880"
GRID = "#e4e3de"

FIGSIZE = (7.0, 4.4)   # inches; ~1400x880 px at 200 dpi
DPI = 200


def apply_style():
    """Install the house rcParams. Call once at the top of a plotting script."""
    mpl.rcParams.update({
        "figure.figsize": FIGSIZE,
        "figure.dpi": 110,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,

        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 12.5,
        "axes.titleweight": "600",
        "axes.titlepad": 30,
        "axes.labelsize": 10.5,
        "axes.labelcolor": INK_2,
        "axes.labelpad": 8,

        # recessive frame: only left + bottom, hairline
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,

        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "grid.alpha": 1.0,

        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 0,
        "ytick.major.size": 0,

        "legend.frameon": False,
        "legend.fontsize": 9.5,
        "legend.handlelength": 1.6,

        "lines.linewidth": 2.0,
        "lines.markersize": 7,
        "lines.markeredgewidth": 1.6,
        "lines.solid_capstyle": "round",
    })


def title(ax, headline, subtitle=None):
    """Left-aligned headline with an optional grey subtitle beneath it."""
    ax.set_title(headline, loc="left", color=INK)
    if subtitle:
        ax.text(0.0, 1.015, subtitle, transform=ax.transAxes,
                fontsize=9.5, color=INK_MUTED, va="bottom", ha="left")


def annotate(ax, x, y, text, dx=6, dy=8, color=INK_2, weight="normal", size=9.0):
    """Direct label next to a data point - used selectively, not on every point."""
    ax.annotate(text, (x, y), textcoords="offset points", xytext=(dx, dy),
                fontsize=size, color=color, fontweight=weight)


def markers(ax, x, y, color, label=None, marker="o"):
    """A line plus ringed markers: a surface-coloured ring keeps overlapping
    marks separable, which matters where two series cross."""
    ax.plot(x, y, color=color, marker=marker, label=label,
            markerfacecolor=color, markeredgecolor=SURFACE, zorder=3)


def save(fig, path_stem):
    """Write both a raster (for the README) and a vector (for the report)."""
    fig.savefig(f"{path_stem}.png")
    fig.savefig(f"{path_stem}.svg")
    plt.close(fig)
    print(f"  wrote {path_stem}.png / .svg")
