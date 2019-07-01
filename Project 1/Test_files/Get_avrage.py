from numpy import load, sum, log10
from os.path import dirname
from sys import path
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from Load_files import get_npz
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_normal, show_plot, save_fig, set_axis_specs
from Quadrature import midpoint, simpson, legendre_gauss
from Test_example import analytical_solution
from matplotlib import pyplot as plt
path.__delitem__(0)


if __name__ == "__main__":
    npzfile = load(rel_path + r'/Data_files/oppgave7d25.npz')
    average= 0
    lambdas = npzfile['lambdas']
    errors_by_lambda = npzfile['errors_by_lambda']
    for i in range(len(lambdas)):
        average+=log10(lambdas[i])*errors_by_lambda[i]
    average = average/sum(errors_by_lambda)
    print ("log10(Average) for d=2.5 =", average, "evt average =", 10**(average))

    npzfile = load(rel_path + r'/Data_files/oppgave7d025.npz')
    average = 0
    lambdas = npzfile['lambdas']
    errors_by_lambda = npzfile['errors_by_lambda']
    for i in range(len(lambdas)):
        average += log10(lambdas[i]) * errors_by_lambda[i]
    average = average / sum(errors_by_lambda)
    print("log10(Average) for d=0.25 =", average, "evt average =", 10**(average))


