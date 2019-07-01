from matplotlib import pyplot as plt
import matplotlib.gridspec as grd
from matplotlib.backends.backend_pdf import PdfPages

def init_plot(navn):
        fig = plt.figure(navn)
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)
        gs = grd.GridSpec(1, 1)

        ax = plt.subplot(gs[0, 0])
        ax.grid(color="lightgrey")

        return fig, ax

def set_axis_specs(ax,x_label,y_label,x_label_size,y_label_size):
    ax.set_xlabel(x_label, fontsize=x_label_size)
    ax.set_ylabel(y_label, fontsize=y_label_size)
    return ax

def plot_logy(ax, x, y, format, label_plot):
    ax.semilogy(x, y, format, label=label_plot, linewidth=2)
    return ax

def plot_logy_test(ax, x, y, format, label_plot):
    ax.semilogy(x, y, format, label=label_plot, linewidth=1.2, markersize=8)
    return ax

def plot_loglog(ax, x, y, format, label_plot):
    ax.loglog(x, y, format, label=label_plot, linewidth = 2)
    return ax

def plot_normal(ax, x, y, format, label_plot):
    ax.plot(x, y, format, label=label_plot, linewidth=2)
    return ax

def plot_small(ax, x, y, format, label_plot):
    ax.plot(x, y, format, label=label_plot, linewidth=1.5)
    return ax

def show_plot(fig, ax):
    ax.legend(fontsize=27)
    plt.show()
    plt.close()

def save_fig(fig, ax, name):
    ax.yaxis.set_label_coords(-0.11, 0.5)
    ax.xaxis.set_label_coords(0.5, -0.07)
    from os.path import dirname
    rel_path = dirname(dirname(__file__))
    with PdfPages(rel_path + "/PDF/"+name) as pdf:
        ax.legend(fontsize=12, loc=1)
        pdf.savefig(fig)

def save_fig_legend_on_top(fig, ax, name):
    ax.yaxis.set_label_coords(-0.10, 0.5)
    ax.xaxis.set_label_coords(0.5, -0.07)
    from os.path import dirname
    rel_path = dirname(dirname(__file__))
    with PdfPages(rel_path + "/PDF/"+name) as pdf:
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=14,
                   bbox_transform=ax.transAxes)
        pdf.savefig(fig)