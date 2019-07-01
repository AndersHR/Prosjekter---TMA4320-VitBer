from numpy import empty, array, linspace, dot, eye, zeros, geomspace
from numpy.linalg import solve
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt

from sys import path
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Quadrature import midpoint, simpson, legendre_gauss
from Load_files import get_npz
path.__delitem__(0)

def fredholm_lhs(d, x_c, x_s, x_q, w_q, f):
    N_c, N_s = x_c.shape[0], x_s.shape[0]
    w_s = get_lagrange_poly_weights(x_s, array([1 for i in range(N_s)]))
    A = empty((N_c, N_s))
    for i in range(N_c):
        for j in range(N_s):
            A[i][j] = num_quadrature(x_q, w_q, get_kernel(f, x_c[i], x_s, w_s, j, d))
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

def get_plot(X, Y, a, b):
    N_plot = 1000
    N_s = len(X)
    w_s = get_lagrange_poly_weights(X, Y)
    x = linspace(a, b, N_plot)
    y = zeros(N_plot)
    for i in range (N_plot):
        for j in range(N_s):
            y[i] += w_s[j]*evaluate_lagrange(X, x[i], j)
    return x, y

if __name__ == "__main__":
    N_q = 900
    x_q, w_q = legendre_gauss(0, 1, N_q)

    fig = plt.figure('Oppgave 2')
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey")
    ax.set_xlabel(r'$x$', fontsize=36)
    ax.set_ylabel(r'$\rho$', fontsize=36)
    c_list_full = ['#ff0000', '#008000', '#0000ff']
    c_list_unsmooth = ['#c30000', '#004900', '#000089']
    c_list_half = ['#f99d9d', '#9ce39c', '#badbff']

    lambdasn = (-4.19044099786)
    lambdas_around = geomspace(10 ** (lambdasn - 0.3), 10 ** (lambdasn + 0.3), 100)
    lambdasn = 10 ** (lambdasn)

    for q_number in range(1, 4):
        a, b, d, x_c, F = get_npz(q_number)
        print ("d = {}".format(d))

        N_s = len(x_c)
        x_s = shift_chebyshev(a, b, chebyshev_n(N_s))
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)


        for i in range(100):
            rho_hat = (solve((dot(A.transpose(), A) + lambdas_around[i] * eye(N_s)), dot(A.transpose(), F)))
            x, y = get_plot(x_s, rho_hat, a, b)
            ax.plot(x, y, c=c_list_half[q_number-1], linestyle='-')

        """rho_hat = (solve((dot(A.transpose(), A) + lambdas_around[99] * eye(N_s)), dot(A.transpose(), F)))
        x, y = get_plot(x_s, rho_hat, a, b)
        ax.plot(x, y, c=c_list_half[q_number-1], linestyle='-', label='Shade q8_{}'.format(q_number))"""
        print('q_number end', q_number)

    for q_number in range(1, 4):
        a, b, d, x_c, F = get_npz(q_number)
        print("d = {}".format(d))

        N_s = len(x_c)
        x_s = shift_chebyshev(a, b, chebyshev_n(N_s))
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)

        rho_hat = (solve((dot(A.transpose(), A) + lambdasn * eye(N_s)), dot(A.transpose(), F)))
        x, y = get_plot(x_s, rho_hat, a, b)

        ax.plot(x, y, 'r-', c=c_list_full[q_number-1], linestyle='-', label=r'$\rho$(x) - q8_{}'.format(q_number), linewidth=4)
        ax.plot(x_s, rho_hat, c=c_list_unsmooth[q_number-1], linestyle='', marker='o', ms=10)

    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=28,
              bbox_transform=ax.transAxes)
    plt.show()
