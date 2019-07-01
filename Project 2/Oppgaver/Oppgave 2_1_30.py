from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, exp, zeros, savez
from numpy.random import seed

def d(T, d_max, s):
    return d_max*exp(-s*T)

if __name__ == "__main__":
    d_max = 30000
    s = 5*1e-4
    num_T = 151
    length = 30

    T_list = linspace(0.00001, 1500, num_T)
    E, E_tot, L_tot, L, last_bend = 0, 0, 0, length, 0
    E_list, L_list = zeros(num_T), zeros(num_T)
    seed_list = [1337, 1, 10352, 9234, 999922, 999999999, 235359, 863286248, 124872, 4652742]

    for k in range (len(seed_list)):

        seed(seed_list[k])
        protein = Protein(length, 0)
        U_mat = protein.potentials

        print("\n{}\n".format(k))

        for i in range (num_T):

            print (T_list[i], d(T_list[i], d_max, s))

            protein = Protein(length, T_list[i])
            protein.potentials = U_mat
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

    file = abs_path + r'/Data_files/Oppgave_2_1_30.npz'
    savez(file, T_list=T_list, E_list=E_list, L_list=L_list)




