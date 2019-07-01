from numpy import empty, array, pi, sin, exp, zeros, geomspace, dot, eye, argmin, savez
from numpy.linalg import solve
from pickle import load
from sys import path
from os.path import dirname
from numpy.random import seed, uniform
from time import time
from matplotlib import pyplot as plt

rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_loglog, show_plot, save_fig, plot_normal, set_axis_specs
from Quadrature import midpoint, simpson, legendre_gauss
from Test_example import analytical_solution
path.__delitem__(0)

def fredholm_rhs(x_c, F):
    N_c = x_c.shape[0]
    b = empty(N_c)
    for i in range(N_c):
        b[i] = F(x_c[i])
    return b

def disturb_rhs(F):
    F_new = F+F*uniform(-1e-3, 1e-3, len(F))
    return F_new


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

def rho(x):
    w, y = 3*pi, -2
    return sin(w*x)*exp(y*x)

if __name__ == "__main__":
    print ('yo')
    d= 0.25
    N_c, N_s = 40, 40
    n_lambda = 200
    N_q = 900



    F_analytic = load(open(rel_path + "/Data_files/F.pkl", "rb"))
    x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
    x_s = shift_chebyshev(0, 1, chebyshev_n(N_s))
    x_q, w_q = legendre_gauss(0, 1, N_q)
    b = F_analytic(x_c, d)
    rho_analytical = rho(x_c)

    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)

    lambdas = geomspace(1e-6, 1.5*1e-3, n_lambda)
    errors_by_lambda = zeros(n_lambda)
    max_error = empty(n_lambda)

    print ("yo^2")
    """for j in range(100000):
        seed(int(time() + 3*j))
        b_tilde = disturb_rhs(b)

        for i in range(n_lambda):
            rho_lambda = solve((dot(A.transpose(), A) + lambdas[i]*eye(N_c)), dot(A.transpose(), b_tilde))
            max_error[i] = max(abs(rho_lambda-rho_analytical))
        errors_by_lambda[argmin(max_error)] += 1

    fig, ax = init_plot('Oppgave 7 -b')
    ax.scatter(lambdas, errors_by_lambda, c='r', s=30)
    ax.set_xscale('log')
    ax.set_xlim(1e-6, 1.5*1e-3)
    ax = set_axis_specs(ax, r'$\lambda$', r'# min |Error|', 23, 23)
    plt.show()

    file = rel_path + r'/Data_files/Oppgave7d025'
    savez(file, lambdas=lambdas, errors_by_lambda=errors_by_lambda)"""

    n_lambda = 100000
    lambdas = geomspace(1e-14, 1, n_lambda)
    max_error = empty(n_lambda)
    seed(int(time()))
    b_tilde = disturb_rhs(b)
    for i in range(n_lambda):
        rho_lambda = solve((dot(A.transpose(), A) + lambdas[i] * eye(N_c)), dot(A.transpose(), b_tilde))
        max_error[i] = max(abs(rho_lambda - rho_analytical))
    fig, ax = init_plot('Oppgave 7 -b')
    plot_loglog(ax, lambdas, max_error, 'r-', "Error($\lambda$)")
    ax = set_axis_specs(ax, r'$\lambda$', r'Error', 30, 30)
    show_plot(fig, ax)











