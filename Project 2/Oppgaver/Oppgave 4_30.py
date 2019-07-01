from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, flip, zeros, append, savez, empty, array, argmin
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt
from numpy.random import seed

if __name__=="__main__":
    seed(2542876)

    bends_per_T = 600
    num_T = 101
    T_min, T_max = 0, 1500
    n_monomers = 30

    T_list = flip(linspace(T_min, T_max, num_T), 0)
    protein = Protein(n_monomers, 0)
    E_list = zeros(1)
    L_list = array([n_monomers])
    bend_list = zeros(1)
    mean_E = empty(num_T)
    mean_L = empty(num_T)

    E, L, last_bend, last_bend_or_change_T = 0, n_monomers, 0, 0

    for T_ind in range(num_T):
        print (T_list[T_ind])
        protein.T=T_list[T_ind]
        E_tot, L_tot = 0, 0
        last_bend_or_change_T = protein.valid_bends
        while (protein.valid_bends < bends_per_T*(T_ind+1)):
            if protein.bend_protein():
                E_list = append(E_list, E)
                L_list = append(L_list, L)
                bend_list = append(bend_list, protein.valid_bends)
                E_tot += E * (protein.valid_bends - last_bend_or_change_T)
                L_tot += L * (protein.valid_bends - last_bend_or_change_T)
                E = protein.u
                L = protein.L
                last_bend = protein.valid_bends
                last_bend_or_change_T = protein.valid_bends
                E_list = append(E_list, E)
                L_list = append(L_list, L)
                bend_list = append(bend_list, protein.valid_bends)
        E_tot += E * (protein.valid_bends - last_bend_or_change_T)
        L_tot += L * (protein.valid_bends - last_bend_or_change_T)
        mean_E[T_ind] = E_tot/bends_per_T
        mean_L[T_ind] = L_tot/bends_per_T
    E_list = append(E_list, E)
    L_list = append(L_list, L)
    bend_list = append(bend_list, protein.valid_bends)

    print (bend_list[argmin(E_list)])

    protein.plot_protein()

    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)

    fig1 = plt.figure('Figur E')
    gs = grd.GridSpec(1, 1)

    ax1 = fig1.add_subplot(gs[0, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.plot(T_list, mean_E, 'r-')
    plt.show()



    file = abs_path + r'/Data_files/Oppgave_4_2.npz'
    savez(file, T_list=T_list, E_list=E_list, L_list=L_list, bend_list=bend_list, mean_E=mean_E, mean_L=mean_L)
