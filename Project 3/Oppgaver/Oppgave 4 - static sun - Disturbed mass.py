from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals
path.__delitem__(0)

from numpy import array, pi, empty, sqrt
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd



G = 4*pi**2

def squared(x, y):
    return sqrt(x**2+y**2)

if __name__=="__main__":
    t_0, t_end = 0, 100
    h = 1e-4
    N_OBJECTS = 2

    m = empty(N_OBJECTS+1)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = get_start_vals('Earth', 1, 0)
    r[1], v[1], m[1] = get_start_vals('Mars', -1, 0)

    m[-1] = 1
    #m[1] = 317.89 / 333480
    #m[0] = 317.89 / 333480
    m[0] = m[-1]

    t, x, y, v_x, v_y = LG4_imp_static_sun(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v, der_v_der_primary, der_v_der_secondary)

    fig1 = plt.figure('Figur r')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)
    ax = fig1.add_subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.set_xlabel(r'$x$   \   Au', fontsize=28)
    ax.set_ylabel(r'$y$   \   Au', fontsize=28)
    ax.set_ylim(-2.5, 2.5)
    ax.set_yticks(array([-1, 0, 1]))
    ax.set_xlim(-2.5, 2.5)
    ax.set_xticks(array([-1, 0, 1]))
    ax.plot(0, 0, 'yo', ms=40, label='Sun')
    ax.plot(x[0], y[0], 'b-', linewidth=1.5, label='Earth')
    ax.plot(x[1], y[1], 'r-', linewidth=1.5, label='Mars')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax.transAxes, markerscale=0.7)

    plt.show()