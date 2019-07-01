from numpy import empty, linspace
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd
from sys import path
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_normal, show_plot, save_fig, set_axis_specs, plot_logy
from Quadrature import midpoint, simpson, legendre_gauss
from Test_example import analytical_solution
path.__delitem__(0)

def F(x, d):
    return ((2/3)-x)/(d*(d**2+(x-2/3)**2)**(1/2)) -((1/3)-x)/(d*(d**2+(x-1/3)**2)**(1/2))

def get_plot(x, d):
    y = empty(len(x))
    for i in range(len(x)):
        y[i] = F(x[i], d)
    return  y

if __name__ == "__main__":
    x = linspace(0, 1, 10**5)
    y_d0025 = get_plot(x, 0.025)
    y_d025 = get_plot(x, 0.25)
    y_d25 = get_plot(x, 2.5)

    fig = plt.figure('Oppgave 2')
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey")
    ax.set_xlabel(r'$x$', fontsize=36)
    ax.set_ylabel("F", fontsize=36)

    ax.semilogy(x, y_d0025, 'r-', label='F(x) - d=0.025', linewidth=4)
    ax.semilogy(x, y_d025, 'b-', label='F(x) - d=0.25', linewidth=4)
    ax.semilogy(x, y_d25, 'g-', label='F(x) - d=2.5', linewidth=4)

    ax.legend(fontsize=28, loc=1)
    plt.show()
    plt.close()