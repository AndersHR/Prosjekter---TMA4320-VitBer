from numpy import empty, array, pi, sin, exp, savez
from numpy.linalg import solve
from pickle import load
from sys import path
from os.path import dirname
from numpy.random import seed, uniform
from matplotlib import pyplot as plt
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from PolyMaker import evaluate_lagrange, get_lagrange_poly_weights
from Chebyshev import chebyshev_n, shift_chebyshev
from Plotting import init_plot, plot_normal, show_plot, save_fig, set_axis_specs, plot_small
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
    N_c, N_s = 40, 40

    seed(1)
    F_analytic = load(open(rel_path + "/Data_files/F.pkl", "rb"))
    x_c = shift_chebyshev(0, 1, chebyshev_n(N_c))
    x_s = shift_chebyshev(0, 1, chebyshev_n(N_s))
    N_q = 900
    x_q, w_q = legendre_gauss(0, 1, N_q)


    #---------------------------------------
    d= 0.025
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    b = F_analytic(x_c, d)
    b_tilde = disturb_rhs(b)

    rho_hat = solve(A, b)
    rho_hat_tilde = solve(A, b_tilde)

    fig1, ax1 = init_plot('Question 6-1 - d=0,025')
    ax1 = set_axis_specs(ax1, r'$x^s$', r'$b$', 18, 18)
    ax1 = plot_normal(ax1, x_c, b, 'r-', r'Unperturbed $b$')
    ax1 = plot_small(ax1, x_c, b_tilde, 'b--', r'Perturbed $\~{b}$')

    fig2, ax2 = init_plot('Question 6-2 - d=0,025')
    ax2 = set_axis_specs(ax2, r'$x^s$', r'$\rho$', 18, 18)
    ax2 = plot_normal(ax2, x_s, rho_hat, 'r-', r'Unperturbed $\hat{\rho}$')
    ax2 = plot_normal(ax2, x_s, rho_hat_tilde, 'b-', r'Perturbed $\hat{\rho}$')
    ax2 = plot_small(ax2, x_s, rho(x_s), 'g--', r'Analytic $\rho$')

    file = rel_path + r'/Data_files/Oppgave6_d0025'
    savez(file, b=b, b_tilde=b_tilde, rho_hat=rho_hat, rho_hat_tilde=rho_hat_tilde)

    #-----------------------------------------------------------------
    d = 0.25
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    b = F_analytic(x_c, d)
    b_tilde = disturb_rhs(b)

    rho_hat = solve(A, b)
    rho_hat_tilde = solve(A, b_tilde)

    fig3, ax3 = init_plot('Question 6-1 - d=0,25')
    ax3 = set_axis_specs(ax3, r'$x^s$', r'$b$', 18, 18)
    ax3 = plot_normal(ax3, x_c, b, 'r-', r'Unperturbed $b$')
    ax3 = plot_small(ax3, x_c, b_tilde, 'b--', r'Perturbed $\~{b}$')

    fig4, ax4 = init_plot('Question 6-2 - d=0,25')
    ax4 = set_axis_specs(ax4, r'$x^s$', r'$\rho$', 18, 18)
    ax4 = plot_normal(ax4, x_s, rho_hat, 'r-', r'Unperturbed $\hat{\rho}$')
    ax4 = plot_normal(ax4, x_s, rho_hat_tilde, 'b-', r'Perturbed $\hat{\rho}$')
    ax4 = plot_small(ax4, x_s, rho(x_s), 'g--', r'Analytic $\rho$')

    file = rel_path + r'/Data_files/Oppgave6_d025'
    savez(file, b=b, b_tilde=b_tilde, rho_hat=rho_hat, rho_hat_tilde=rho_hat_tilde)

    #------------------------------------------------------------------------
    d = 2.5
    A = fredholm_lhs(d, x_c, x_s, x_q, w_q, f)
    b = F_analytic(x_c, d)
    b_tilde = disturb_rhs(b)

    rho_hat = solve(A, b)
    rho_hat_tilde = solve(A, b_tilde)

    fig5, ax5 = init_plot('Question 6-1 - d=2,5')
    ax5.set_xlim(0, 0.2)
    ax5 = set_axis_specs(ax5, r'$x^s$', r'$b$', 18, 18)
    ax5 = plot_normal(ax5, x_c, b, 'r-', r'Unperturbed $b$')
    ax5 = plot_small(ax5, x_c, b_tilde, 'b--', r'Perturbed $\~{b}$')

    fig6, ax6 = init_plot('Question 6-2 - d=2,5')
    ax6 = set_axis_specs(ax6, r'$x^s$', r'$\rho$', 18, 18)
    ax6 = plot_normal(ax6, x_s, rho_hat, 'r-', r'Unperturbed $\hat{\rho}$')
    ax6 = plot_normal(ax6, x_s, rho_hat_tilde, 'b-', r'Perturbed $\hat{\rho}$')
    ax6 = plot_small(ax6, x_s, rho(x_s), 'g--', r'Analytic $\rho$')

    file = rel_path + r'/Data_files/Oppgave6_d25'
    savez(file, b=b, b_tilde=b_tilde, rho_hat=rho_hat, rho_hat_tilde=rho_hat_tilde)


    ax1.legend(fontsize=20)
    ax2.legend(fontsize=20)
    ax3.legend(fontsize=20)
    ax4.legend(fontsize=20)
    ax5.legend(fontsize=20)
    ax6.legend(fontsize=20)

    plt.show()





