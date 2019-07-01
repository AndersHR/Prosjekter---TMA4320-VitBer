from numpy import empty, array, pi, sin, exp, linspace, zeros
from numpy.linalg import solve
from pickle import load
from sys import path
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_normal, show_plot, save_fig
from Quadrature import midpoint, simpson, legendre_gauss
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

def get_biggest_difference(x_s, rho_hat, x_list):
    w_s = get_lagrange_poly_weights(x_s, rho_hat)
    ylist = zeros(len(x_list))
    for i in range(len(x_list)):
        for j in range(len(x_s)):
            ylist[i] += w_s[j]*evaluate_lagrange(x_s, x_list[i], j)

    return ylist





if __name__ == "__main__":
    d = 0.025
    n_min, n_max = 25, 40
    size = 100
    N_q = 300


    N_c_list = array([i for i in range(n_min, n_max+1)])
    x_q, w_q = legendre_gauss(0, 1, N_q)
    F_analytic = load(open(rel_path + "/Data_files/F.pkl", "rb"))
    max_diff, max_diff_rho = empty(len(N_c_list)), empty(len(N_c_list))
    y_to_plot = empty((n_max-n_min+1, size))
    x_list = linspace(0, 1, size)

    for i in N_c_list:
        print (i)
        N_c, N_s = i, i
        x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
        x_s = x_c
        A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
        b = F_analytic(x_c, d)
        rho_hat = solve(A, b)

        y_to_plot[i-n_min] = get_biggest_difference(x_s, rho_hat, x_list)
        max_diff_rho[i-n_min] = max(abs(y_to_plot[i-n_min]-rho(x_list)))
        max_diff[i-n_min] = max(abs(rho_hat-rho(x_s)))
    fig, ax = init_plot('Question 5')
    plot_normal(ax, N_c_list, max_diff_rho, 'b-', r'Err whole intervall')
    plot_normal(ax, N_c_list, max_diff, 'r-', r'Err(N$_c$)')
    show_plot(fig, ax)




