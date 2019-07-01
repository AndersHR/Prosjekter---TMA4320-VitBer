from numpy import empty, array, linspace, meshgrid
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as grd
from matplotlib import cm


from sys import path
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_normal, show_plot, save_fig, plot_logy_test, save_fig_legend_on_top, set_axis_specs, plot_small
from Quadrature import midpoint, simpson, legendre_gauss
from scipy.linalg import solve
from Test_example import analytical_solution

path.__delitem__(0)


def fredholm_lhs(d, x_c, x_s, x_q, w_q, f):
    N_c, N_s = x_c.shape[0], x_s.shape[0]
    w_s = get_lagrange_poly_weights(x_s, array([1 for i in range(N_s)]))
    A = empty((N_c, N_s))
    for i in range(N_c):
        for j in range(N_s):
            A[i][j] = num_quadrature(x_q, w_q, get_kernel(f, x_c[i], x_s, w_s, j, d))
            #A[i][j] = num_quadrature(x_q, w_q, x_c[i], x_s, d, j)
    return A

def f(x, y, d):
    return d/(d**2+(y-x)**2)**(3/2)

def get_kernel(f, x, x_s, w_s, j, d):
    return lambda y: f(x, y, d)*w_s[j]*evaluate_lagrange(x_s, y, j)

def num_quadrature(x_q, w_q, K):
    sum_q = 0
    n_q = len(x_q)
    for i in range(n_q):
        sum_q += K(x_q[i])*w_q[i]
    return sum_q

if __name__ == "__main__":
    N_q = 900
    x_q, w_q = legendre_gauss(0, 1, N_q)
    print ("hei1")
    d = 0.25
    N_c, N_s = 30, 30
    x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
    x_s = x_c
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    print (A[0]-A[4])
    fig = plt.figure('3d')
    plt.rc('xtick', labelsize=17)
    plt.rc('ytick', labelsize=17)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0], projection='3d')
    ax.set_xlabel('j', fontsize=26)
    ax.set_ylabel('i', fontsize=36)
    x, y = linspace(1, 30, 30), linspace(1, 30, 30)
    X, Y = meshgrid(x, y)

    colors = cm.get_cmap('Reds')(A)
    rcount, ccount, _ = colors.shape

    surf = ax.plot_wireframe(X, Y, A, color='r', linewidth=3)
    surf.set_facecolor((0, 0, 0, 0))
    plt.show()


