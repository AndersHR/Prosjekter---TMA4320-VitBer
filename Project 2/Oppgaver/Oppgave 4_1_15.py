from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, flip, zeros, append, savez, array
from numpy.random import seed

if __name__=="__main__":
    seed(3567238576)

    bends_per_T = 600
    num_T = 51
    T_min, T_max = 0.0000000000001, 1500
    n_monomers = 15

    T_list = flip(linspace(T_min, T_max, num_T), 0)
    protein = Protein(n_monomers, 0)
    E_list = zeros(1)
    L_list = array([n_monomers])
    bend_list = zeros(1)


    E, L, last_bend = 0, n_monomers, 0

    for T_ind in range(num_T):
        print (T_list[T_ind])
        protein.T=T_list[T_ind]
        E_tot, L_tot = 0, 0
        while (protein.valid_bends < bends_per_T*(T_ind+1)):
            if protein.bend_protein():
                E_list = append(E_list, E)
                L_list = append(L_list, L)
                bend_list = append(bend_list, protein.valid_bends)
                E = protein.u
                L = protein.L
                last_bend = protein.valid_bends
                E_list = append(E_list, E)
                L_list = append(L_list, L)
                bend_list = append(bend_list, protein.valid_bends)
    E_list = append(E_list, E)
    L_list = append(L_list, L)
    bend_list = append(bend_list, protein.valid_bends)

    file_data = abs_path + r'/Data_files/Oppgave_4_1_15.npz'
    savez(file_data, T_list=T_list, E_list=E_list, L_list=L_list, bend_list=bend_list, protein=protein.aa_pos,
          protein_vals=[protein.length, protein.u, protein.L])

    file_pot = abs_path + r'/Data_files/Oppgave_4_1_15_pot.npz'
    savez(file_pot, potentials=protein.potentials)