from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import exp, zeros, append
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt
from numpy.random import seed

def d(T, d_max, s):
    return d_max*exp(-s*T)


if __name__ == "__main__":
    seed(76567)
    d_max = 10000
    s = 1e-4

    E, E_tot_0, E_tot_5, last_bend = 0, 0, 0, 0
    E_list_0, bend_list_0 = zeros(1), zeros(1)
    E_list_5, bend_list_5 = zeros(1), zeros(1)

    protein_0 = Protein(15, 1e-15)
    protein_5 = Protein(15, 500)
    protein_5.potentials = protein_0.potentials

    it =0

    while (protein_0.valid_bends < 5000):
        it+=1
        print (protein_0.valid_bends, it)
        if protein_0.bend_protein():
            E_list_0 = append(E_list_0, E)
            E_tot_0 += E * (protein_0.valid_bends - last_bend)
            last_bend = protein_0.valid_bends
            E = protein_0.u
            E_list_0 = append(E_list_0, E)
            bend_list_0 = append(bend_list_0, [last_bend, last_bend])
    E_tot_0 += E * (protein_0.valid_bends - last_bend)
    E_tot_0 = E_tot_0/protein_0.valid_bends
    E_list_0 = append(E_list_0, E)
    bend_list_0 = append(bend_list_0, protein_0.valid_bends)

    E, last_bend = 0, 0

    while (protein_5.valid_bends < 5000):
        print (protein_5.valid_bends)
        if protein_5.bend_protein():
            E_list_5 = append(E_list_5, E)
            E_tot_5 += E * (protein_5.valid_bends - last_bend)
            last_bend = protein_5.valid_bends
            E = protein_5.u
            E_list_5 = append(E_list_5, E)
            bend_list_5 = append(bend_list_5, [last_bend, last_bend])
    E_tot_5 += E * (protein_5.valid_bends - last_bend)
    E_tot_5 = E_tot_5/protein_0.valid_bends
    E_list_5 = append(E_list_5, E)
    bend_list_5 = append(bend_list_5, protein_5.valid_bends)


    fig = plt.figure('Figur E')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(100, 1)

    ax1 = plt.subplot(gs[:49, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_xlim(-1, 5001)
    ax1.set_ylim(-6.2 * 1e-20, 0.1 * 1e-20)
    ax1.set_yticks([-6e-20, -4e-20, -2e-20, 0])
    ax1.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    ax1.set_ylabel (r'$E$    /    J', fontsize = 28)
    ax1.yaxis.set_label_coords(-0.03, -0.05)
    k100_plot, = ax1.plot(bend_list_0, E_list_0, 'r-', label = r'$T=10^{-15}$K')
    ax1.tick_params(
        axis = 'x',
        which='both',
        labelbottom = 'off'
    )

    ax2 = plt.subplot(gs[50:, 0])
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.set_xlim(-1, 5001)
    ax2.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    ax2.set_ylim(-7.05 * 1e-20, 0.2 * 1e-20)
    ax2.set_yticks([-6e-20, -4e-20, -2e-20, 0])
    ax2.yaxis.get_offset_text().set_visible(False)
    ax2.set_xlabel(r'Number of valid bends', fontsize = 28)
    k500_plot, = ax2.plot(bend_list_5, E_list_5, 'b-', label=r'$T=500$K')

    ax1.legend(handles = [k100_plot, k500_plot], bbox_to_anchor=(0.1, 1.02, 0.9, .102), loc=3, ncol=2, mode="expand", borderaxespad=0., fontsize=30,
               bbox_transform=ax1.transAxes)

    plt.show()




