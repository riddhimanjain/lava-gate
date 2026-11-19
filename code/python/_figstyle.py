"""_figstyle.py -- shared house style for the paper's figures."""
import matplotlib as mpl
import matplotlib.pyplot as plt

C_BLUE   = "#1f77b4"
C_ORANGE = "#ff7f0e"
C_GREEN  = "#2ca02c"
C_RED    = "#d62728"
C_PURPLE = "#9467bd"
C_GREY   = "#7f7f7f"

BLACK    = "black"
BAND     = "#eeeeee"

CHROM_A, CHROM_B = "#8fb8d6", C_BLUE

C_POS, C_NEG = C_RED, C_BLUE

MS_SERIES   = 4
S_SCATTER   = 22
A_SCATTER   = 0.75
A_GRID      = 0.3
FS_LEGEND   = 8

LW_REF      = 1.0
LW_ZERO     = 0.8
LW_EMPH     = 2.0
A_QUIET     = 0.6

S_FILLED    = 22
S_HOLLOW    = 20
MEW         = 0.9
S_DENSE     = 5
LW_INTERVAL = 0.9
A_INTERVAL  = 0.45
LW_RUG      = 0.4
A_RUG       = 0.75
A_BAND      = 0.18

PANEL_W, PANEL_H = 5.0, 4.4


def apply():
    """Install the template's look. Call once at the top of a figure script."""
    mpl.rcParams.update({
        "font.family":        "sans-serif",
        "font.sans-serif":    ["DejaVu Sans"],
        "font.size":          10.0,
        "mathtext.fontset":   "dejavusans",

        "axes.titlesize":     12.0,
        "axes.labelsize":     10.0,
        "xtick.labelsize":    10.0,
        "ytick.labelsize":    10.0,
        "legend.fontsize":    FS_LEGEND,

        "axes.spines.top":    True,
        "axes.spines.right":  True,
        "axes.linewidth":     0.8,
        "axes.edgecolor":     "black",

        "axes.grid":          True,
        "axes.axisbelow":     True,
        "grid.alpha":         A_GRID,
        "grid.linestyle":     "-",
        "grid.linewidth":     0.8,

        "legend.frameon":     True,
        "legend.framealpha":  0.8,

        "lines.linewidth":    1.5,
        "lines.markersize":   MS_SERIES,

        "figure.facecolor":   "white",
        "axes.facecolor":     "white",
        "savefig.facecolor":  "white",

        "figure.dpi":         110,
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
    })


def refline(ax, value, axis="y", kind="identity", zorder=2):
    """A reference line, drawn as the template draws them."""
    if kind == "identity":
        kw = dict(color=BLACK, ls="--", lw=LW_REF)
    else:
        kw = dict(color=BLACK, lw=LW_ZERO)
    kw["zorder"] = zorder
    return (ax.axhline if axis == "y" else ax.axvline)(value, **kw)


def save(fig, path):
    """PNG (as the template) plus vector PDF for a journal."""
    fig.tight_layout()
    fig.savefig(f"{path}.png")
    fig.savefig(f"{path}.pdf")
    plt.close(fig)
    print(f"  wrote {path.split('/')[-1]}.png + .pdf")
