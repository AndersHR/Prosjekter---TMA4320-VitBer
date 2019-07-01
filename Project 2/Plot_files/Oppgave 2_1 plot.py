from os.path import dirname

abs_path = dirname(dirname(__file__))
from numpy import load, zeros
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt

if __name__=="__main__":
    T_list, E_list, L_list = zeros(151), zeros(151), zeros(151)

    file1 = abs_path + r'/Data_files/Oppgave_2_1_15_test.npz'
    npz_file1 = load(file1)
    T_list += npz_file1['T_list']
    E_list += npz_file1['E_list']
    L_list += npz_file1['L_list']
    npz_file1.close()

    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)

    fig1 = plt.figure('Figur E')
    gs = grd.GridSpec(1, 1)

    ax1 = fig1.add_subplot(gs[0, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_xlabel(r'$T$    /    K', fontsize=24)
    ax1.set_ylabel(r'$\langle E \rangle$   /    J', fontsize=28)
    ax1.set_xlim(-30, 1530)
    ax1.set_xticks([0, 250, 500, 750, 1000, 1250, 1500])
    ax1.set_ylim(-6.1e-20, -1.8e-20)
    ax1.set_yticks([-6e-20, -5e-20, -4e-20, -3e-20, -2e-20,])
    ax1.plot(T_list, E_list, '--', c='#FF7575', linewidth=1.5)
    ax1.plot(T_list, E_list, 'o', c='#FF0000', ms=6)
    plt.show()

    fig2 = plt.figure('Figur E')
    gs = grd.GridSpec(1, 1)

    ax2 = fig2.add_subplot(gs[0, 0])
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.set_xlabel(r'$T$    /    K', fontsize=24)
    ax2.set_ylabel(r'$\langle L \rangle$   /    Number of monomers', fontsize=28)
    ax2.set_xlim(-30, 1530)
    ax2.set_xticks([0, 250, 500, 750, 1000, 1250, 1500])
    ax2.plot(T_list, L_list, '--', c='#FF7575', linewidth=1.5)
    ax2.plot(T_list, L_list, 'o', c='#FF0000', ms=6)
    plt.show()
