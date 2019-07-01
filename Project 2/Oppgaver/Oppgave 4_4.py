from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import savez, empty, load, copy
from numpy.random import seed

if __name__=="__main__":
    seed(6443578)

    T = 510
    n_monomers = 30
    max_bends = 1e6

    npzfile = load(abs_path + r'/Data_files/Oppgave_4_1_{}_pot.npz'.format(str(n_monomers)))
    potentials = copy(npzfile['potentials'])
    npzfile.close()

    protein = Protein(n_monomers, T)
    lowest_E = 0
    lowest_L = 0
    protein_lowest_E = empty(n_monomers)

    protein.potentials = potentials
    print (potentials-protein.potentials)

    E, L, last_bend_or_change_T = 0, n_monomers, 0

    while (protein.valid_bends < max_bends):
        if (protein.valid_bends%1000==0):
            print (protein.valid_bends)
        if protein.bend_protein():
            if (protein.u < lowest_E):
                lowest_E = protein.u
                protein_lowest_E = copy(protein.aa_pos)
                print (protein_lowest_E)
                lowest_L = protein.L



    file = abs_path + r'/Data_files/Oppgave_4_4_{}.npz'.format(str(n_monomers))
    savez(file, protein_vals=[n_monomers, lowest_E, lowest_L], protein = protein_lowest_E)
