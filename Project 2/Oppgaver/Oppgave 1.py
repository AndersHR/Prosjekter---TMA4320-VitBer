from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy.random import seed

if __name__ == "__main__":
    seed (27643)

    protein = Protein(10,300)

    protein.plot_protein()

    protein.bend_protein()

    protein.plot_protein()

    while(protein.valid_bends < 2):
        protein.bend_protein()

    protein.plot_protein()


