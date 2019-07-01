from numpy import load, flip
from os.path import dirname
rel_path = dirname(dirname(__file__))

def get_npz(q_number):
    f = open(rel_path + r'\Data_files\q8_' + str(q_number) + '.npz', 'rb')
    npzfile = load(f)
    a = npzfile['a']
    b = npzfile['b']
    d = npzfile['d']
    x_c = flip(npzfile['xc'], 0)
    F = flip(npzfile['F'], 0)
    f.close()
    return a, b, d, x_c, F

