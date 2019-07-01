from numpy import empty, array, pi, sin, exp, dot, arange, savez, inf
from numpy.linalg import norm
from pickle import load
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

def fredholm_rhs(x_c, F):
    N_c = x_c.shape[0]
    b = empty(N_c)
    for i in range(N_c):
        b[i] = F(x_c[i])
    return b

def fredholm_lhs(d, x_c, x_s, w_s, x_q, w_q, f):
    N_c, N_s = x_c.shape[0], x_s.shape[0]
    A = empty((N_c, N_s))
    for i in range(N_c):
        for j in range(N_s):
            A[i][j] = num_quadrature(x_q, w_q, get_integrand(f, x_c[i], x_s, w_s, j, d))
    return A

def f(x, y, d):
    return d/(d**2+(y-x)**2)**(3/2)

def get_integrand(f, x, x_s, w_s, j, d):
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
    d=0.025
    n_c = 40
    n_s = 40

    x_c = shift_chebyshev(0, 1, chebyshev_n(n_c))
    x_s = shift_chebyshev(0, 1, chebyshev_n(n_s))
    w_s = get_lagrange_poly_weights(x_s, array([1 for i in range(n_s)]))
    rho_hat = array([rho(x_s[i]) for i in range (n_s)])
    F_analytic = load(open(rel_path + "/Data_files/F.pkl", "rb"))
    f_eval = F_analytic(x_c, d)

    fig, ax = init_plot('Oppgave 3')
    ax = set_axis_specs(ax,'$N_q$','Max Error',24,24)

    panel_list = arange(20, 301, 1)
    y_list_midpoint = empty(len(panel_list))
    y_list_simpson = empty(len(panel_list))
    y_list_LG = empty(len(panel_list))
    for i in range(len(panel_list)):
        x_q, w_q = midpoint(0, 1, panel_list[i])
        A = fredholm_lhs(d, x_c, x_s, w_s, x_q, w_q, f)
        F = dot(A, rho_hat)
        y_list_midpoint[i] = norm((F - f_eval), 0)

    ax = plot_logy(ax, panel_list, y_list_midpoint, 'r-', "Err($N_q$) - midpoint")

    for i in range(len(panel_list)):
        x_q, w_q = simpson(0, 1, panel_list[i])
        A = fredholm_lhs(d, x_c, x_s, w_s, x_q, w_q, f)
        F = dot(A, rho_hat)
        y_list_simpson[i] = norm((F - f_eval), inf)

    ax = plot_logy(ax, panel_list, y_list_simpson, 'g-', "Err($N_q$) - Simpson")

    for i in range(len(panel_list)):
        x_q, w_q = legendre_gauss(0, 1, panel_list[i])
        A = fredholm_lhs(d, x_c, x_s, w_s, x_q, w_q, f)
        F = dot(A, rho_hat)
        y_list_LG[i] = norm(F - f_eval, inf)

    ax = plot_logy(ax, panel_list, y_list_LG, 'b-', "Err($N_q$) - L_G")

    show_plot(fig, ax)

    file = rel_path + r'/Data_files/Oppgave2'
    savez(file, panel_list=panel_list, midpoint=y_list_midpoint, simpson=y_list_simpson, LG=y_list_LG)

