from numpy import load
from sys import path
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from Plotting import init_plot, plot_normal, show_plot, save_fig, plot_logy_test, save_fig_legend_on_top, set_axis_specs, plot_small
path.__delitem__(0)

if __name__ == "__main__":
    npzfile = load(rel_path+r'/Data_files/Oppgave2.npz')
    fig = plt.figure('Oppgave 2')
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey")
    ax.set_xlabel(r'$N_q$', fontsize=36)
    ax.set_ylabel("Max Error", fontsize=36)

    ax.semilogy(npzfile['panel_list'], npzfile['midpoint'], 'r-', label=r'Err($N_q$) - Midpoint', linewidth=4)
    ax.semilogy(npzfile['panel_list'], npzfile['simpson'], 'g-', label=r'Err($N_q$) - Simpson', linewidth=4)
    ax.semilogy(npzfile['panel_list'], npzfile['LG'], 'b-', label=r'Err($N_q$) - Legendre-Gauss', linewidth=4)

    ax.legend(bbox_to_anchor=(0.605, 0.89, 0.39, 0.102),loc=1, mode="expand", borderaxespad=0., fontsize=26,
                   bbox_transform=ax.transAxes)
    plt.show()
    plt.close()