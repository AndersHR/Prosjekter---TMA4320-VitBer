from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))


from numpy import load, copy, array
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd

if __name__=="__main__":
    cont_file = abs_path + r'/Data_files/Oppgave_3_2_cont_v2.npz'
    cont_npzfile = load(cont_file)
    delta_vs_cont, sum_delta_vs_cont = copy(cont_npzfile['delta_vs']), copy(cont_npzfile['sum_delta_vs'])
    cont_npzfile.close()

    Hohmann_file = abs_path + r'/Data_files/Oppgave_3_2_Hohmann.npz'
    Hohmann_npzfile = load(Hohmann_file)
    n_transfers_Hohmann, sum_delta_vs_Hohmann = copy(Hohmann_npzfile['n_transfers'])[0], copy(Hohmann_npzfile['sum_delta_vs'])
    Hohmann_npzfile.close()

    fig = plt.figure('Continous acceleration - semilogx')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel(r'$a$   \   $\frac{\mathrm{m}}{\mathrm{s}}$', fontsize=28)
    ax.set_ylabel(r'$\sum{\Delta v}$   \   $\frac{\mathrm{km}}{\mathrm{s}}$', fontsize=28)
    ax.set_ylim(4, 17)
    ax.set_yticks(array([5, 7, 9, 11, 13, 15]))
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.semilogx(delta_vs_cont, sum_delta_vs_cont*1e-3, 'ro--', label = "Continous acceleration", linewidth=2, ms=10)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand",
              borderaxespad=0., fontsize=22,
              bbox_transform=ax.transAxes)
    plt.show()

    fig = plt.figure('Hohmann-transfers')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel(r'Number of Hohmann-transfers', fontsize=28)
    ax.set_ylabel(r'$\sum{\Delta v}$   \   $\frac{\mathrm{km}}{\mathrm{s}}$', fontsize=28)
    ax.set_ylim(3.7, 4.7)
    ax.set_yticks(array([3.8, 4.0, 4.2, 4.4, 4.6]))
    ax.set_xticks(array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.plot(array([i for i in range(1, n_transfers_Hohmann+1)]), sum_delta_vs_Hohmann * 1e-3, 'ro--', label="Hohmann transfers", linewidth=2, ms=10)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand",
              borderaxespad=0., fontsize=22,
              bbox_transform=ax.transAxes)
    plt.show()

    print (delta_vs_cont[0], sum_delta_vs_cont[0])
    print (sum_delta_vs_Hohmann[-1])
