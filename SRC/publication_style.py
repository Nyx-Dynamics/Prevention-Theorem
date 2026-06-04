"""
Publication-quality matplotlib style helpers.

Provides set_publication_style(), despine(), and get_colorblind_palette()
used by the stochastic PEP analysis scripts.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt


def set_publication_style():
    """Apply publication-quality defaults to matplotlib."""
    mpl.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'lines.linewidth': 1.5,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
    })


def despine(ax, top=True, right=True, left=False, bottom=False):
    """Remove specified spines from an axes (seaborn-style)."""
    if top:
        ax.spines['top'].set_visible(False)
    if right:
        ax.spines['right'].set_visible(False)
    if left:
        ax.spines['left'].set_visible(False)
    if bottom:
        ax.spines['bottom'].set_visible(False)


def get_colorblind_palette(n=8):
    """Return a colorblind-safe palette of n colors (Wong 2011 + extensions)."""
    base = [
        '#000000',  # black
        '#E69F00',  # orange
        '#56B4E9',  # sky blue
        '#009E73',  # green
        '#F0E442',  # yellow
        '#0072B2',  # blue
        '#D55E00',  # vermillion
        '#CC79A7',  # purple-pink
    ]
    palette = (base * ((n // len(base)) + 1))[:n]
    return palette
