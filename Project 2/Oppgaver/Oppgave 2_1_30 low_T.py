from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import load, exp, zeros, savez, append, copy
from numpy.random import seed, randint
from time import time

def d(T, d_max, s):
    return d_max*exp(-s*T)

if __name__ == "__main__":
    d_max = 30000
    s = 5*1e-4
    num_T = 151
    length = 30

    file1 = abs_path + r'/Data_files/Oppgave_2_1_30_test.npz'
    npz_file1 = load(file1)

    T_list, E_list, L_list = zeros(151), zeros(151), zeros(151)
    T_list += npz_file1['T_list']
    E_list += npz_file1['E_list']
    L_list += npz_file1['L_list']
    seed_list = copy(npz_file1['seed_list'])

    npz_file1.close()

    for i in range(51):
        E_list[i] = E_list[i]*(len(seed_list))
        L_list[i] = L_list[i]*(len(seed_list))

    E, E_tot, L_tot, L, last_bend = 0, 0, 0, length, 0

    start = time()

    while (len(seed_list) < 51):

        new_seed = randint(10, int(4294967290/2))*2-randint(0, 2)
        valid_seed = True
        for old_seed in seed_list:
            if old_seed == new_seed:
                valid_seed = False

        if valid_seed:
            print ("\n{}\n".format(len(seed_list)+1))
            seed((new_seed))
            seed_list = append(seed_list, new_seed)
            protein0 = Protein(length, 0)

            for i in range (51):

                print (T_list[i], d(T_list[i], d_max, s))

                protein = Protein(length, T_list[i])
                protein.potentials = protein0.potentials
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

    print (time()-start)

    for i in range(51):
        E_list[i] = E_list[i]/(len(seed_list))
        L_list[i] = L_list[i]/(len(seed_list))


    file = abs_path + r'/Data_files/Oppgave_2_1_30_test.npz'
    savez(file, T_list=T_list, E_list=E_list, L_list=L_list, seed_list=seed_list)




