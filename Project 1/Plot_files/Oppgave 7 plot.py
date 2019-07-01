from numpy import load, pi, sin, exp, linspace
from sys import path
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from Plotting import init_plot, plot_normal, show_plot, plot_loglog,save_fig, plot_logy_test, save_fig_legend_on_top, set_axis_specs, plot_small
path.__delitem__(0)

if __name__ == "__main__":
    npzfile_d025 = load(rel_path + r'/Data_files/Oppgave7d025.npz')
    npzfile_d25 = load(rel_path + r'/Data_files/Oppgave7d25.npz')

    fig, ax = init_plot('Oppgave 7 -b')
    ax.scatter(npzfile_d25['lambdas'][1:], npzfile_d25['errors_by_lambda'][1:], c='r', s=30)
    ax.set_xscale('log')
    ax.set_xlim(1e-11, 3 * 1e-9)
    ax = set_axis_specs(ax, r'$\lambda$', r'# min |Error|', 30, 30)
    plt.show()



    fig, ax = init_plot('Oppgave 7 -b')
    ax.scatter(npzfile_d025['lambdas'][1:], npzfile_d025['errors_by_lambda'][1:], c='r', s=30)
    ax.set_xscale('log')
    ax.set_xlim(1e-6, 1.5 * 1e-3)
    ax = set_axis_specs(ax, r'$\lambda$', r'# min |Error|', 30, 30)
    plt.show()