from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, flip, savez, empty, load, copy
from numpy.random import seed

if __name__=="__main__":
    seed(3567238576)

    bends_per_T = 3000
    num_T = 151
    T_min, T_max = 0.0000000000001, 1500
    n_monomers = 15

    npzfile = load(abs_path + r'/Data_files/Oppgave_4_1_15_pot.npz')
    potentials = copy(npzfile['potentials'])
    npzfile.close()

    T_list = flip(linspace(T_min, T_max, num_T), 0)
    protein = Protein(n_monomers, 0)
    mean_E = empty(num_T)
    mean_L = empty(num_T)

    print (potentials-protein.potentials)

    E, L, last_bend_or_change_T = 0, n_monomers, 0

    for T_ind in range(num_T):
        print (T_list[T_ind])
        protein.T=T_list[T_ind]
        E_tot, L_tot = 0, 0
        last_bend_or_change_T = protein.valid_bends
        while (protein.valid_bends < bends_per_T*(T_ind+1)):
            if protein.bend_protein():
                E_tot += E * (protein.valid_bends - last_bend_or_change_T)
                L_tot += L * (protein.valid_bends - last_bend_or_change_T)
                E = protein.u
                L = protein.L
                last_bend_or_change_T = protein.valid_bends
        E_tot += E * (protein.valid_bends - last_bend_or_change_T)
        L_tot += L * (protein.valid_bends - last_bend_or_change_T)
        mean_E[T_ind] = E_tot/bends_per_T
        mean_L[T_ind] = L_tot/bends_per_T




    file = abs_path + r'/Data_files/Oppgave_4_2_15_151T_3000bends.npz'
    savez(file, T_list=T_list, mean_E=mean_E, mean_L=mean_L)
