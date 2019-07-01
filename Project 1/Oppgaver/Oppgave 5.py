from numpy import empty, array, pi, sin, exp, savez
import numpy.linalg
from matplotlib import pyplot as plt
from pickle import load
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

def fredholm_rhs(x_c, F):
    N_c = x_c.shape[0]
    b = empty(N_c)
    for i in range(N_c):
        b[i] = F(x_c[i])
    return b

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
    n_min, n_max = 5, 30

    N_c_list = array([i for i in range(n_min, n_max+1)])
    d_list = array([0.025, 0.25, 2.5])
    c_list = array(['r-', 'b-', 'g-'])

    N_q = 900
    x_q, w_q = legendre_gauss(0, 1, N_q)

    max_diff_1, max_diff_2, max_diff_3 = empty(len(N_c_list)), empty(len(N_c_list)), empty(len(N_c_list))

    F_analytic = load(open(rel_path + "/Data_files/F.pkl", "rb"))

    fig, ax = init_plot('Question 5 - Error d=0.025 & 0.25')
    d=0.025
    for i in N_c_list:
        print (i)
        N_c, N_s = i, i
        x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
        x_s = x_c
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
        b = F_analytic(x_c, d)
        rho_hat = solve(A, b)
        max_diff_1[i-n_min] = max(abs(rho_hat-rho(x_c)))

    plot_logy_test(ax, N_c_list, max_diff_1, 'ro--', r'Err(N$_{c}$) - d = ' + '{}'.format(d))

    d = 0.25
    for i in N_c_list:
        print(i)
        N_c, N_s = i, i
        x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
        x_s = x_c
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
        b = F_analytic(x_c, d)
        rho_hat = numpy.linalg.solve(A, b)


        max_diff_2[i - n_min] = max(abs(rho_hat - rho(x_c)))

    plot_logy_test(ax, N_c_list, max_diff_2, 'bo--', r'Err($N_{c}$) - d = ' + '{}'.format(d))
    ax = set_axis_specs(ax, r'$N_c$', r'Max $|$Error$_j|$', 24, 24)
    #save_fig_legend_on_top(fig, ax, 'Question 5 - Error d=0,025 & 0,25.pdf')

    #fig, ax = init_plot('Question 5 - Error d=2.5')
    d = 2.5
    for i in N_c_list:
        print(i)
        N_c, N_s = i, i
        x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
        x_s = x_c
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
        b = F_analytic(x_c, d)
        rho_hat = numpy.linalg.solve(A, b)
        max_diff_3[i - n_min] = max(abs(rho_hat - rho(x_c)))

    plot_logy_test(ax, N_c_list, max_diff_3, 'go--', r'Err($N_c$) - d ={}'.format(d))
    show_plot(fig, ax)
    #ax = set_axis_specs(ax, r'$N_c$', r'Max $|$Error$_j|$', 18, 18)
    #save_fig_legend_on_top(fig, ax, 'Question 5 - Error d=2,5.pdf')
    file = rel_path + r'/Data_files/Oppgave5'
    savez(file, N_c_list=N_c_list, max_diff_1=max_diff_1, max_diff_2=max_diff_2, max_diff_3=max_diff_3)


    d = 0.25
    N_c, N_s = 29,29
    x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
    x_s = x_c
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    b = F_analytic(x_c, d)
    rho_hat = numpy.linalg.solve(A, b)

    fig, ax = init_plot('Oppgave 5 - rho_hat for d=0.25 & N_c= 30')

    ax = plot_normal(ax, x_s, rho_hat, 'ro--', r'$\^{\rho}$')
    ax = plot_small(ax, x_s, rho(x_s), 'b-', r'Analytic $\rho$')

    ax = set_axis_specs(ax, r'$x^s$', r'$\rho$', 18, 18)

    ax.text(0.95, 0.95, r'$N_c$ = {}'.format(N_c), fontsize=15, horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=0.5))

    save_fig_legend_on_top(fig, ax, 'Question 5-1-2.pdf')
    #save_fig(fig, ax, 'Question 5 - d=0,25 & N_c = 29.pdf')

    plt.show()
    plt.close()

    d = 2.5
    N_c, N_s = 9, 9
    x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
    x_s = x_c
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    b = F_analytic(x_c, d)
    rho_hat = numpy.linalg.solve(A, b)

    fig, ax = init_plot('Oppgave 5 - rho_hat for d=0.25 & N_c= 30')
    ax = plot_normal(ax, x_s, rho_hat, 'ro--', r'$\^{\rho}$')
    ax = plot_small(ax, x_s, rho(x_s), 'b-', r'Analytic $\rho$')
    ax = set_axis_specs(ax, r'$x^s$', r'$\rho$', 18, 18)

    ax.text(0.95, 0.95, r'$N_c$ = {}'.format(N_c), fontsize=15, horizontalalignment='right', verticalalignment='top',
             transform=ax.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))

    save_fig_legend_on_top(fig, ax, 'Question 5-2-2.pdf')
    #save_fig(fig, ax, 'Question 5 - d=2.5 & N_c = 9.pdf') """

    plt.show()
    plt.close()






