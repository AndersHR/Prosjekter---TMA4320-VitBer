from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, exp, zeros, savez
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt
from numpy.random import seed

def d(T, d_max, s):
    return d_max*exp(-s*T)

if __name__ == "__main__":
    d_max = 20000
    s = 5*1e-4
    num_T = 151
    length = 15

    T_list = linspace(0.00001, 1500, num_T)
    E, E_tot, L_tot, L, last_bend = 0, 0, 0, length, 0
    E_list, L_list = zeros(num_T), zeros(num_T)
    old_seed_list = [1337, 1, 10352, 9234, 999922, 999999999, 235359, 863286248, 124872, 4652742]
    seed_list = [186072, 895242, 418811, 153930, 446805, 511978, 537066, 436767, 457513, 739423]
    for k in range (len(seed_list)):
        print ("\n{}\n".format(k))
        seed(seed_list[k])
        for i in range (num_T):
            print (T_list[i], d(T_list[i], d_max, s))
            protein = Protein(length, T_list[i])
            E, E_tot, L_tot, L, last_bend = 0, 0, 0, length, 0
            while (protein.valid_bends < d(T_list[i], d_max, s)):
                if protein.bend_protein():
                    L_tot += L*(protein.valid_bends-last_bend)
                    E_tot += E*(protein.valid_bends-last_bend)
                    last_bend = protein.valid_bends
                    L = protein.L
                    E = protein.u
            E_tot += E * (protein.valid_bends - last_bend)
            L_tot += L * (protein.valid_bends - last_bend)
            E_list[i] += E_tot/protein.valid_bends
            L_list[i] += L_tot/protein.valid_bends
    for i in range(num_T):
        E_list[i] = E_list[i]/len(seed_list)
        L_list[i] = L_list[i]/len(seed_list)

    """fig = plt.figure('Figur E')
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.plot(T_list, E_list, 'r-')
    plt.show()

    ig = plt.figure('Figur E')
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.plot(T_list, L_list, 'r-')
    plt.show()"""

    file = abs_path + r'/Data_files/Oppgave_2_1_15_v2.npz'
    savez(file, T_list=T_list, E_list=E_list, L_list=L_list)




