from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import array, exp, savez, float64, amin, dot, empty
from numpy.random import uniform, seed
from scipy.constants import k
from time import time

def d(T, d_max, s):
    return d_max*exp(-s*T)


if __name__ == "__main__":
    arr1 = array([0.001,1e-50, 0])
    arr2 = empty(3)
    T=0
    for i in range(3):
        arr2[i] = exp(-arr1[i]/T)
    print (arr2)
    """protein15 = Protein(15, 500)
    protein30 = Protein(30, 500)
    start = time()
    it = 0
    while (protein15.valid_bends<100000):
        protein15.bend_protein()
        it+=1
    time15 = time()-start
    bends_tried15=it
    start = time()
    it = 0
    while (protein30.valid_bends < 100000):
        protein30.bend_protein()
        it += 1
    time30 = time() - start
    bends_tried30=it
    print ("Tid 15: {}\nBends_15: {}\nBends tried 15: {}\n\n".format(time15, protein15.valid_bends, bends_tried15))
    print ("Tid 30: {}\nBends_30: {}\nBends tried 30: {}\n\n".format(time30, protein30.valid_bends, bends_tried30))"""

